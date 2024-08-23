#!/usr/bin/env python3
from argparse import ArgumentParser
from itertools import chain
from os import walk
from pathlib import Path
from typing import Iterable, Mapping, Set

import barcodeutils as bu
import manhole
from common import BARCODE_UMI_FASTQ_PATH, TRANSCRIPT_FASTQ_PATH
from fastq_utils import Read, fastq_reader, find_grouped_fastq_files, revcomp

BARCODE_LENGTH = 16
BARCODE_STARTS = [0]
BARCODE_SEGMENTS = [slice(start, start + BARCODE_LENGTH) for start in BARCODE_STARTS]
UMI_SEGMENT = slice(16, 28)

BARCODE_QUAL_DUMMY = "F" * BARCODE_LENGTH * len(BARCODE_STARTS)


class KeyDefaultDict(dict):
    def __missing__(self, key):
        return key


def read_barcode_allowlist(barcode_filename: Path) -> Set[str]:
    print("Reading barcode allowlist from", barcode_filename)
    with open(barcode_filename) as f:
        return set(f.read().split())


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


def main(
    metadata_dir: Path,
    fastq_dirs: Iterable[Path],
    output_dir: Path = Path(),
):
    output_dir.mkdir(exist_ok=True, parents=True)

    visium_version_number = get_visium_plate_version(metadata_dir)

    barcode_filename = f"/opt/visium-v{visium_version_number}.txt"

    barcode_allowlist = read_barcode_allowlist(barcode_filename)
    correcter = bu.BarcodeCorrecter(barcode_allowlist, edit_distance=1)

    buf = output_dir / BARCODE_UMI_FASTQ_PATH
    trf = output_dir / TRANSCRIPT_FASTQ_PATH

    all_fastqs = chain.from_iterable(
        find_grouped_fastq_files(fastq_dir, 2) for fastq_dir in fastq_dirs
    )

    with open(buf, "w") as cbo, open(trf, "w") as tro:
        for barcode_umi_fastq, transcript_fastq in all_fastqs:
            usable_count = 0
            i = 0
            print("Correcting barcodes in", transcript_fastq, "and", barcode_umi_fastq)
            transcript_reader = fastq_reader(transcript_fastq)
            barcode_umi_reader = fastq_reader(barcode_umi_fastq)
            for i, (tr, br) in enumerate(zip(transcript_reader, barcode_umi_reader), 1):
                barcode_pieces = [br.seq[s] for s in BARCODE_SEGMENTS]
                corrected = [correcter.correct(barcode) for barcode in barcode_pieces]
                if all(corrected):
                    usable_count += 1
                    umi_seq = br.seq[UMI_SEGMENT]
                    umi_qual = br.qual[UMI_SEGMENT]
                    new_seq = "".join(corrected + [umi_seq])
                    new_qual = BARCODE_QUAL_DUMMY + umi_qual
                    new_br = Read(
                        read_id=br.read_id,
                        seq=new_seq,
                        unused=br.unused,
                        qual=new_qual,
                    )
                    print(tr.serialize(), file=tro)
                    print(new_br.serialize(), file=cbo)

            print("Total count:", i)
            print("Usable count:", usable_count)
            print("Proportion:", usable_count / i)


if __name__ == "__main__":
    manhole.install(activate_on="USR1")

    p = ArgumentParser()
    p.add_argument("metadata_dir", type=Path)
    p.add_argument("fastq_dirs", type=Path, nargs="+")
    args = p.parse_args()

    main(
        metadata_dir=args.metadata_dir,
        fastq_dirs=args.fastq_dirs,
    )
