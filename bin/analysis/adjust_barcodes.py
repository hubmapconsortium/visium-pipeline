#!/usr/bin/env pypy3
from argparse import ArgumentParser
from pathlib import Path
from typing import Iterable

import correct_visium_barcodes
import manhole
from common import ADJ_OUTPUT_DIR, Assay


def main(metadata_dir: Path, input_dirs: Iterable[Path]):
    correct_visium_barcodes.main(metadata_dir, input_dirs, output_dir=ADJ_OUTPUT_DIR)


if __name__ == "__main__":
    manhole.install(activate_on="USR1")

    p = ArgumentParser()
    p.add_argument("metadata_dir", type=Path)
    p.add_argument("directory", type=Path, nargs="+")
    args = p.parse_args()

    main(args.metadata_dir, args.directory)
