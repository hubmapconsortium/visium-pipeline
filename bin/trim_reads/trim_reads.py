#!/usr/bin/env python3
from argparse import ArgumentParser
from concurrent.futures import ProcessPoolExecutor, wait
from itertools import chain
from os import fspath
from pathlib import Path
from shlex import quote
from shutil import copy
from subprocess import check_call
from typing import Iterable, Sequence, Tuple

from fastq_utils import find_grouped_fastq_files, fastq_reader, Read

from common import (
    BARCODE_UMI_FASTQ_PATH,
    TRANSCRIPT_FASTQ_GZ_PATH,
    TRANSCRIPT_FASTQ_PATH,
    Assay,
)

OUTPUT_PATH = Path("trimmed")

TRIM_COMMAND = [
    "/opt/seqtk",
    "trimfq",
    "{input_fastq}",
]


def find_adj_fastq_files(directory: Path) -> Tuple[Path, Path]:
    # not general enough to implement in fastq-utils; very specific
    # to how we create "synthetic" barcode + UMI FASTQ files
    barcode_umi_fastq = directory / BARCODE_UMI_FASTQ_PATH

    transcript_fastq = directory / TRANSCRIPT_FASTQ_PATH
    transcript_fastq_gz = directory / TRANSCRIPT_FASTQ_GZ_PATH

    if transcript_fastq.is_file():
        return barcode_umi_fastq, transcript_fastq
    elif transcript_fastq_gz.is_file():
        return barcode_umi_fastq, transcript_fastq_gz
    else:
        message = (
            f"Couldn't find {TRANSCRIPT_FASTQ_PATH} or {TRANSCRIPT_FASTQ_GZ_PATH} in {directory}"
        )
        raise ValueError(message)


def trim_reads(fastq_r1: Path, fastq_r2: Path, output_subdir: Path):
    print("Copying", quote(fspath(fastq_r1)), "to", quote(fspath(output_subdir)))
    copy(fastq_r1, output_subdir)

    command = [piece.format(input_fastq=fastq_r2) for piece in TRIM_COMMAND]
    fastq_r2_out = output_subdir / fastq_r2.name
    command_str = " ".join(quote(s) for s in command)
    print("Running", command_str, "with output", quote(fspath(fastq_r2_out)))
    with open(fastq_r2_out, "wb") as f:
        check_call(command, stdout=f)

def trim_seq_and_qual(read):
    return Read(read.read_id, read.seq[0:50], read.unused, read.qual[0:50])

def trim_reads_visium_ffpe(fastq_r1: Path, fastq_r2: Path, output_subdir: Path):
    print("Copying", quote(fspath(fastq_r1)), "to", quote(fspath(output_subdir)))
    copy(fastq_r1, output_subdir)

    reads = fastq_reader(fastq_r2)
    trimmed_reads = [trim_seq_and_qual(read) for read in reads]
    fastq_r2_out = output_subdir / fastq_r2.name
    with open(fastq_r2_out, "a") as f:
        for read in trimmed_reads:
            f.write(read.serialize())
            f.write('\n')

def main(adj_fastq_dir: Path, threads: int):
    fastq_pairs: Iterable[Sequence[Path]]
    fastq_pairs = [find_adj_fastq_files(adj_fastq_dir)]

    with ProcessPoolExecutor(max_workers=threads) as executor:
        futures = []
        for i, (r1_fastq_file, r2_fastq_file) in enumerate(fastq_pairs, 1):
            subdir = OUTPUT_PATH / str(i)
            subdir.mkdir(exist_ok=True, parents=True)
            callable = trim_reads_visium_ffpe
            future = executor.submit(
                callable,
                r1_fastq_file,
                r2_fastq_file,
                subdir,
            )
            futures.append(future)
        wait(futures)


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("adj_fastq_dir", type=Path)
    p.add_argument("-p", "--threads", type=int)
    args = p.parse_args()

    main(args.adj_fastq_dir, args.threads)
