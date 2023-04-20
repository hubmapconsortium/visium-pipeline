#!/usr/bin/env python3
import csv
import re
from argparse import ArgumentParser
from os import environ, fspath, walk
from pathlib import Path
from subprocess import check_call, check_output, CalledProcessError
from typing import Iterable, Optional, Sequence, Tuple
import pandas as pd

import simplesam
import manhole
import anndata
from fastq_utils import find_grouped_fastq_files

from common import (
    BARCODE_UMI_FASTQ_PATH,
    TRANSCRIPT_FASTQ_GZ_PATH,
    TRANSCRIPT_FASTQ_PATH,
    Assay,
)

base_index = '/opt/gencode.v35.intron-exon.sidx'
base_transcript_map = '/opt/gencode.v35.annotation.expanded.tx2gene.tsv'

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

def get_visium_probe_set_version(directory: Path, probe_set_version_parameter: int = None)-> int:
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
                probe_set_version_metadata = int(metadata[metadata_probe_set_version_field])
                print(
                    f"Read expected cell count from {maybe_metadata_file}: {probe_set_version_metadata}"
                )

    present_probe_set_versions = sum(x is not None for x in [probe_set_version_parameter, probe_set_version_metadata])
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

def read_expected_cell_count(directory: Path) -> Optional[int]:
    cell_count_from_file = None
    cell_count_metadata = None

    cell_count_file = directory / cell_count_filename
    if cell_count_file.is_file():
        with open(cell_count_file) as f:
            cell_count_from_file = int(f.read().strip())
            print(f"Read expected cell count from {cell_count_file}: {cell_count_from_file}")

    maybe_metadata_file = find_metadata_file(directory)
    if maybe_metadata_file and maybe_metadata_file.is_file():
        with open(maybe_metadata_file, newline="") as f:
            r = csv.DictReader(f, delimiter="\t")
            metadata = next(r)
            if (
                metadata_cell_count_field in metadata
                and metadata[metadata_cell_count_field].isdigit()
            ):
                cell_count_metadata = int(metadata[metadata_cell_count_field])
                print(
                    f"Read expected cell count from {maybe_metadata_file}: {cell_count_metadata}"
                )

    present_cell_counts = sum(x is not None for x in [cell_count_from_file, cell_count_metadata])
    if present_cell_counts == 0:
        return None
    elif present_cell_counts == 1:
        return cell_count_from_file or cell_count_metadata or 0
    else:
        if cell_count_from_file == cell_count_metadata:
            return cell_count_from_file
        else:
            message = (
                f"Found mismatched cell counts: {cell_count_from_file} in {cell_count_file}, "
                f"and {cell_count_metadata} in {maybe_metadata_file}"
            )
            raise ValueError(message)


def read_expected_cell_counts(directories: Sequence[Path]) -> Optional[int]:
    cell_counts = []
    for directory in directories:
        cell_count = read_expected_cell_count(directory)
        if cell_count is not None:
            cell_counts.append(cell_count)

    dirs_with_cell_counts = len(cell_counts)
    if dirs_with_cell_counts == 0:
        return None
    elif dirs_with_cell_counts == len(directories):
        total_expected = sum(cell_counts)
        print("Total expected cells:", total_expected)
        return total_expected
    else:
        message = (
            f"Found expected cell counts in {dirs_with_cell_counts} of "
            f"{len(directories)} directories, need 0 or {len(directories)} "
            f"input directories with cell counts (can't mix auto-detection "
            f"and guided cell barcode counting)"
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
    assay: Assay,
    orig_fastq_dirs: Sequence[Path],
    trimmed_fastq_dir: Path,
    expected_cell_count: Optional[int],
    keep_all_barcodes: bool,
    threads: Optional[int],
    visium_probe_set_version: Optional[int],
):
    threads = threads or 1

    visium_plate_version = get_visium_plate_version(orig_fastq_dirs[0])
    visium_probe_set_version = get_visium_probe_set_version(orig_fastq_dirs[0], visium_probe_set_version)

    index = f"/opt/v{visium_probe_set_version}.fasta"

    fastq_pairs = list(find_grouped_fastq_files(trimmed_fastq_dir, 2))

    if not fastq_pairs:
        raise ValueError("No FASTQ files found")

    # hack
    if assay in {Assay.VISIUM_FF, Assay.VISIUM_FFPE}:
        # Don't support multiple input directories for Slide-seq; this will
        # likely cause significantly incorrect results due to barcode overlap
        # between multiple input data sets
        if len(orig_fastq_dirs) != 1:
            raise ValueError("Need exactly 1 input directory for Visium")

    if assay in {Assay.VISIUM_FFPE, Assay.VISIUM_FF}:
        barcode_file = f'/opt/visium-v{visium_plate_version}.txt'
        r1_fastq_file, r2_fastq_file = fastq_pairs[0]
        BWA_INDEX_COMMAND = f"bwa index {index}"
        UMI_EXTRACT_COMMAND = f"umi_tools extract --extract-method=string --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNNNN --stdin {r1_fastq_file} --stdout extracted_barcode_umi.fastq.gz --read2-in={r2_fastq_file} --read2-out=extracted_transcript.tar.gz"
        check_call(UMI_EXTRACT_COMMAND, shell=True)

        check_call(BWA_INDEX_COMMAND)
        BWA_COMMAND = f"bwa mem -M -t {threads} {index} {r1_fastq_file} {r2_fastq_file} > out.sam"
        check_call(BWA_COMMAND)

        with simplesam.Reader(open('out.sam')) as in_sam:
            with simplesam.Writer(open('mapped.sam', 'w')) as out_sam:
                for read in in_sam:
                    if read.mapped:
                        ensembl_id = str(read).split('\t')[2]
                        read['XT'] = ensembl_id
                        out_sam.write(read)

        SAMTOOLS_COMMAND = f"samtools view -S -b mapped.sam > mapped.bam && samtools sort -@ {threads} mapped.bam -o sorted.bam && samtools index sorted.bam"
        check_call(SAMTOOLS_COMMAND)
        UMI_DEDUP_COMMAND = "umi_tools count --per-gene --gene-tag=XT --per-cell -I assigned_sorted.bam -S counts.tsv.gz"
        check_call(UMI_DEDUP_COMMAND)

        adata = anndata.read_umi_tools("counts.tsv.gz")
        adata.write("expr.h5ad")

if __name__ == "__main__":
    manhole.install(activate_on="USR1")

    p = ArgumentParser()
    p.add_argument("assay", choices=list(Assay), type=Assay)
    p.add_argument("trimmed_fastq_dir", type=Path)
    p.add_argument("orig_fastq_dir", type=Path, nargs="+")
    p.add_argument("--expected-cell-count", type=int)
    p.add_argument("--keep-all-barcodes", action="store_true")
    p.add_argument("-p", "--threads", type=int)
    p.add_argument("--visium-probe-set-version", type=int, nargs="?")
    args = p.parse_args()

    main(
        args.assay,
        args.orig_fastq_dir,
        args.trimmed_fastq_dir,
        args.expected_cell_count,
        args.keep_all_barcodes,
        args.threads,
        args.visium_probe_set_version,
    )
