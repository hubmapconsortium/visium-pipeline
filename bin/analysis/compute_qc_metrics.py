#!/usr/bin/env python3
import json
import logging
import HTSeq
from argparse import ArgumentParser
from pathlib import Path
from typing import List, Optional, Tuple
from collections import Counter
from os import fspath

import numpy as np
import anndata
import manhole
import pandas as pd
import scanpy as sc

from common import Assay


def write_scanpy_qc(adata: anndata.AnnData):
    qc_by_cell, qc_by_gene = sc.pp.calculate_qc_metrics(adata)

    qc_path = Path("qc_results.hdf5").absolute()
    print("Saving QC results to", qc_path)
    with pd.HDFStore(qc_path) as store:
        store["qc_by_cell"] = qc_by_cell
        store["qc_by_gene"] = qc_by_gene


def write_alignment_qc(bam_file):

    bam = HTSeq.BAM_Reader(fspath(bam_file))

    barcode_counts = Counter()
    barcode_reads_in_peaks = Counter()

    total_reads = 0
    mapped_reads = 0
    unmapped_reads = 0
    alignment_qualities = []

    i = 0
    for i, seg in enumerate(bam):
        if not (i % 10000):
            logging.debug(f"Processed {i} reads")
        barcode = seg.read.name.split(":", 1)[1]
        barcode_counts[barcode] += 1

        total_reads += 1
        if seg.aligned:
            mapped_reads += 1
            alignment_qualities.append(seg.aQual)
        else:
            unmapped_reads += 1
            continue

    logging.info(f"Processed peak overlap for {i} reads")


    proportion_mapped = mapped_reads / total_reads

    percentiles = [
        (0, "minimum"),
        (25, "25th_quantile"),
        (50, "median"),
        (75, "75th_quantile"),
        (100, "maximum"),
    ]
    five_number_dict = {key: np.percentile(alignment_qualities, pc) for pc, key in percentiles}

    qc_report = {
        "total_reads": total_reads,
        "mapped_reads": mapped_reads,
        "unmapped_reads": unmapped_reads,
        "mapped_proportion": proportion_mapped,
        "mapping_quality": five_number_dict,
    }

    output_file = Path("alignment_qc_report.json")
    logging.info(f"Writing QC measures to {output_file}")
    with open(output_file, "w") as text_file:
        json.dump(qc_report, text_file, indent=4)

def main(assay: Assay, h5ad_primary: Path, bam_file: Path):
    expr_primary = anndata.read_h5ad(h5ad_primary)
    if assay.secondary_analysis_layer in expr_primary.layers:
        expr_primary.X = expr_primary.layers[assay.secondary_analysis_layer]
    expr_primary.var_names_make_unique()

    write_scanpy_qc(expr_primary)
    write_alignment_qc(bam_file)


if __name__ == "__main__":
    manhole.install(activate_on="USR1")
    p = ArgumentParser()
    p.add_argument("assay", choices=list(Assay), type=Assay)
    p.add_argument("h5ad_primary", type=Path)
    p.add_argument("bam_file", type=Path)
    args = p.parse_args()

    main(args.assay, args.h5ad_primary, args.bam_file)
