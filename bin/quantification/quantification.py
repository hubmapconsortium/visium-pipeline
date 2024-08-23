#!/usr/bin/env python3
import csv
import re
from argparse import ArgumentParser
from os import environ, fspath, walk
from pathlib import Path
from subprocess import CalledProcessError, check_call, check_output
from typing import Iterable, Optional, Sequence, Tuple

import anndata
import manhole
import pandas as pd
import simplesam
from common import (
    BARCODE_UMI_FASTQ_PATH,
    TRANSCRIPT_FASTQ_GZ_PATH,
    TRANSCRIPT_FASTQ_PATH,
    Assay,
)
from fastq_utils import find_grouped_fastq_files

base_index = "/opt/gencode.v35.intron-exon.sidx"
base_transcript_map = "/opt/gencode.v35.annotation.expanded.tx2gene.tsv"

cell_count_filename = "extras/expected_cell_count.txt"
metadata_filename_pattern = re.compile(r"^[0-9A-Fa-f]{32}-metadata.tsv$")
metadata_cell_count_field = "expected_cell_count"
metadata_probe_set_version_field = "visium_probe_set_version"
barcode_whitelist_path = Path("barcode_whitelist.txt")


def find_metadata_file(directory: Path) -> Optional[Path]:
    """
    Finds and returns the first metadata file for a HuBMAP data set.
    Does not check whether the dataset ID (32 hex characters) matches
    the directory name, nor whether there might be multiple metadata files.
    """
    for file_path in directory.iterdir():
        if metadata_filename_pattern.match(file_path.name):
            return file_path


def find_files(directory: Path, pattern: str) -> Iterable[Path]:
    for dirpath_str, dirnames, filenames in walk(directory):
        dirpath = Path(dirpath_str)
        for filename in filenames:
            filepath = dirpath / filename
            if filepath.match(pattern):
                yield filepath


def get_visium_plate_version(directory: Path) -> int:
    gpr_file = list(find_files(directory, "*.gpr"))[0]
    return int(gpr_file.stem[1])


def get_visium_probe_set_version(
    directory: Path, probe_set_version_parameter: int = None
) -> int:
    probe_set_version_metadata = None
    maybe_metadata_file = find_metadata_file(directory)
    if maybe_metadata_file and maybe_metadata_file.is_file():
        with open(maybe_metadata_file, newline="") as f:
            r = csv.DictReader(f, delimiter="\t")
            metadata = next(r)
            if (
                metadata_probe_set_version_field in metadata
                and metadata[metadata_probe_set_version_field].isdigit()
            ):
                probe_set_version_metadata = int(
                    metadata[metadata_probe_set_version_field]
                )
                print(
                    f"Read expected cell count from {maybe_metadata_file}: {probe_set_version_metadata}"
                )

    present_probe_set_versions = sum(
        x is not None for x in [probe_set_version_parameter, probe_set_version_metadata]
    )
    if present_probe_set_versions == 0:
        return None
    elif present_probe_set_versions == 1:
        return probe_set_version_parameter or probe_set_version_metadata or 0
    else:
        if probe_set_version_parameter == probe_set_version_metadata:
            return probe_set_version_parameter
        else:
            message = (
                f"Found mismatched probe set versions: {probe_set_version_parameter} as parameter, "
                f"and {probe_set_version_metadata} in {maybe_metadata_file}"
            )
            raise ValueError(message)


def find_adj_fastq_files(directory: Path) -> Iterable[Tuple[Path, Path]]:
    # not general enough to implement in fastq-utils; very specific
    # to how we create "synthetic" barcode + UMI FASTQ files
    for subdir in directory.iterdir():
        barcode_umi_fastq = subdir / BARCODE_UMI_FASTQ_PATH

        transcript_fastq = subdir / TRANSCRIPT_FASTQ_PATH
        transcript_fastq_gz = subdir / TRANSCRIPT_FASTQ_GZ_PATH

        if transcript_fastq.is_file():
            yield barcode_umi_fastq, transcript_fastq
        elif transcript_fastq_gz.is_file():
            yield barcode_umi_fastq, transcript_fastq_gz


def main(
    trimmed_fastq_dir: Path,
    metadata_dir: Path,
    expected_cell_count: Optional[int],
    keep_all_barcodes: bool,
    threads: Optional[int],
    visium_probe_set_version: Optional[int],
):
    threads = threads or 1

    visium_plate_version = get_visium_plate_version(metadata_dir)
    visium_probe_set_version = get_visium_probe_set_version(
        metadata_dir, visium_probe_set_version
    )

    index = f"/opt/v{visium_probe_set_version}.fasta"
    copy_command = f"ln -s {index} v{visium_probe_set_version}.fasta"
    check_call(copy_command, shell=True)
    index = f"v{visium_probe_set_version}.fasta"

    fastq_pairs = list(find_adj_fastq_files(trimmed_fastq_dir))

    if not fastq_pairs:
        raise ValueError("No FASTQ files found")

    barcode_file = f"/opt/visium-v{visium_plate_version}.txt"
    r1_fastq_file, r2_fastq_file = fastq_pairs[0]
    BWA_INDEX_COMMAND = f"bwa index {index}"
    check_call(BWA_INDEX_COMMAND, shell=True)
    UMI_EXTRACT_COMMAND = f"umi_tools extract --extract-method=string --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNNNN --stdin {r1_fastq_file} --stdout extracted_barcode_umi.fastq.gz --read2-in={r2_fastq_file} --read2-out=extracted_transcript.fastq.gz"
    check_call(UMI_EXTRACT_COMMAND, shell=True)

    BWA_COMMAND = (
        f"bwa mem -M -t {threads} {index} extracted_transcript.fastq.gz > out.sam"
    )
    check_call(BWA_COMMAND, shell=True)

    SAMTOOLS_COMMAND = f"samtools view -S -b -t /opt/v{visium_probe_set_version}.fasta.fai out.sam > out.bam && samtools sort -@ {threads} out.bam -o sorted.bam && samtools index sorted.bam"

    check_call(SAMTOOLS_COMMAND, shell=True)

    UMI_DEDUP_COMMAND = (
        "umi_tools count --per-contig --per-cell -I sorted.bam -S counts.tsv.gz"
    )
    check_call(UMI_DEDUP_COMMAND, shell=True)

    adata = anndata.read_umi_tools("counts.tsv.gz")
    adata.write("expr.h5ad")


if __name__ == "__main__":
    manhole.install(activate_on="USR1")

    p = ArgumentParser()
    p.add_argument("trimmed_fastq_dir", type=Path)
    p.add_argument("metadata_dir", type=Path)
    p.add_argument("--expected-cell-count", type=int)
    p.add_argument("--keep-all-barcodes", action="store_true")
    p.add_argument("-p", "--threads", type=int)
    p.add_argument("--visium-probe-set-version", type=int, nargs="?")
    args = p.parse_args()

    main(
        args.trimmed_fastq_dir,
        args.metadata_dir,
        args.expected_cell_count,
        args.keep_all_barcodes,
        args.threads,
        args.visium_probe_set_version,
    )
