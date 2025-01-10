#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path
from subprocess import check_call
from typing import Tuple

import manhole
from fastq_utils import collect_fastq_files_by_directory

READ_REMOVAL_COMMAND_TEMPLATE = [
    "/opt/scrubber/scripts/scrub.sh",
    "-i",
    "{input_path}",
    "-o",
    "{out_path}",
    "-p",
    "{threads}"
]


def unzip_input_file(input_file: Path, unzipped_dir) -> Path:
    if ".gz" in input_file.suffixes:
        unzip_command = f"gunzip -c {input_file} > {unzipped_dir}/{input_file.name}"
        print(unzip_command)
        check_call(unzip_command, shell=True)
        return unzipped_dir / Path(input_file.stem)
    else:
        return input_file

def zip_output_file(output_file: Path):
    zip_command = f"gzip {output_file}"
    check_call(zip_command, shell=True)

def single_file_human_read_remove(fastq_file_and_subdirs: Tuple[Path, Path, Path], threads: int):
    """
    Run human read removal on a single fastq file

    Takes an absolute path to the input file and a relative path to
    the output subdirectory
    """

    input_path = fastq_file_and_subdirs[0]
    unzipped_input_path = unzip_input_file(input_path, fastq_file_and_subdirs[2])
    output_path = fastq_file_and_subdirs[1] / fastq_file_and_subdirs[0].name
    command = [piece.format(input_path=unzipped_input_path, out_path=output_path, threads=threads) for piece in READ_REMOVAL_COMMAND_TEMPLATE]
    print("Running", " ".join(command))
    check_call(command, shell=True)
    zip_output_file(output_path)

def main(directory: Path, threads: int):
    """
    Crawl directory, create appropriate output subdirectories based on input directory structure
    """
    fastq_files_by_directory = collect_fastq_files_by_directory(directory)
    print("Found", len(fastq_files_by_directory), "directories containing FASTQ files")

    out_dir = Path("sanitized_fastqs")
    unzipped_dir = Path("unzipped")

    for directory, files in fastq_files_by_directory.items():
        subdir = out_dir / directory
        unzipped_subdir = unzipped_dir / directory
        subdir.mkdir(exist_ok=True, parents=True)
        unzipped_subdir.mkdir(exist_ok=True, parents=True)

        for file in files:
            single_file_human_read_remove((file, subdir, unzipped_subdir), threads)


if __name__ == "__main__":
    manhole.install(activate_on="USR1")

    p = ArgumentParser()
    p.add_argument("directory", type=Path)
    p.add_argument("threads", type=int)
    args = p.parse_args()

    main(args.directory, args.threads)
