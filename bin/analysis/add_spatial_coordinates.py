#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path
from typing import Optional

import anndata
import manhole
import pandas as pd
from common import Assay
from os import walk
from typing import Iterable
import read_visium_positions
import numpy as np

barcode_matching_dir = "barcode_matching"

def apply_affine_transform(coords, affine):
    ones = np.ones((coords.shape[0],1))
    coords_concat = np.append(coords, ones, axis=1)
    coords_transform = coords_concat @ affine
    coords_trim = np.delete(coords_transform, 2, axis=1)
    return coords_trim

def find_files(directory: Path, pattern: str) -> Iterable[Path]:
    for dirpath_str, dirnames, filenames in walk(directory):
        dirpath = Path(dirpath_str)
        for filename in filenames:
            filepath = dirpath / filename
            if filepath.match(pattern):
                yield filepath

def annotate(h5ad_path: Path, dataset_dir: Path, assay: Assay, img_dir: Optional[Path] = None, metadata_dir: Optional[Path] = None) -> anndata.AnnData:
    barcode_pos, slide_id, scale_factor, spot_diameter, affine_matrix = read_visium_positions.read_visium_positions(metadata_dir, img_dir)
    d.obs['Tissue Coverage Fraction'] = barcode_pos['Tissue Coverage Fraction']
#        d.obs['Tissue'] = barcode_pos['Tissue']
    spatial_key = "spatial"
    library_id = slide_id
    d.uns[spatial_key] = {library_id: {}}
    d.uns[spatial_key][library_id]["scalefactors"] = {"tissue_hires_scalef": 1.0, "spot_diameter_fullres": spot_diameter}

    quant_bc_set = set(d.obs.index)
    pos_bc_set = set(barcode_pos.index)
    overlap = quant_bc_set & pos_bc_set
    positions_overlap = barcode_pos.loc[list(overlap), :]

    quant_minus_pos = quant_bc_set - pos_bc_set
    positions_missing = pd.DataFrame(
        index=list(quant_minus_pos),
        columns=barcode_pos.columns,
        dtype=float,
    )

    quant_pos = pd.concat([positions_overlap, positions_missing])
    quant_pos_ordered = quant_pos.loc[d.obs.index, ["X", 'Y']] if assay in {assay.VISIUM_FF} else quant_pos.loc[d.obs.index]

    d.obsm["X_spatial"] = apply_affine_transform(quant_pos_ordered.to_numpy(), affine_matrix)
    d.obsm["spatial"] = d.obsm["X_spatial"]
    d.obsm["X_spatial_gpr"] = quant_pos_ordered.to_numpy()


if __name__ == "__main__":
    manhole.install(activate_on="USR1")

    p = ArgumentParser()
    p.add_argument("h5ad_path", type=Path)
    p.add_argument("base_dir", type=Path)
    p.add_argument("assay", choices=list(Assay), type=Assay)
    args = p.parse_args()

    d = annotate(args.h5ad_path, args.base_dir, args.assay)
    d.write_h5ad("spatial_pos_annotated.h5ad")
