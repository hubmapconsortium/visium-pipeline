#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path
from typing import Optional, Sequence

import anndata
import manhole

import add_spatial_coordinates
from common import Assay

H5AD_PATH = Path("expr.h5ad")


def main(
    assay: Assay,
    h5ad_file: Path,
    raw_fastq_dirs: Sequence[Path],
    img_dir: Optional[Path],
    metadata_dir: Optional[Path],
    metadata_json: Optional[Path],
):

    expr_data = add_spatial_coordinates.annotate(h5ad_file, raw_fastq_dirs[0], assay, img_dir, metadata_dir)
    print(expr_data.obs.columns)
    expr_data.write_h5ad(H5AD_PATH)


if __name__ == "__main__":
    manhole.install(activate_on="USR1")

    p = ArgumentParser()
    p.add_argument("assay", choices=list(Assay), type=Assay)
    p.add_argument("h5ad_file", type=Path)
    p.add_argument("raw_fastq_dir", type=Path, nargs="+")
    p.add_argument("--img_dir", type=Path)
    p.add_argument("--metadata_dir", type=Path)
    p.add_argument("--metadata_json", type=Path)
    args = p.parse_args()

    print(args)

    main(args.assay, args.h5ad_file, args.raw_fastq_dir, args.img_dir, args.metadata_dir, args.metadata_json)
