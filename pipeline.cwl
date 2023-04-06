#!/usr/bin/env cwl-runner

class: Workflow
cwlVersion: v1.0
label: scRNA-seq pipeline using Salmon and Alevin
requirements:
  SubworkflowFeatureRequirement: {}
  ScatterFeatureRequirement: {}
inputs:
  fastq_dir:
    label: "Directory containing FASTQ files"
    type: Directory[]
  assay:
    label: "scRNA-seq assay"
    type: string
  threads:
    label: "Number of threads for Salmon"
    type: int
    default: 1
  expected_cell_count:
    type: int?
  keep_all_barcodes:
    type: boolean?
  visium_probe_set_version:
    type: int?
outputs:
  count_matrix_h5ad:
    outputSource: salmon_quantification/count_matrix_h5ad
    type: File
    label: "Unfiltered count matrix from Alevin, converted to H5AD, spliced and unspliced counts"
  fastqc_dir:
    outputSource: fastqc/fastqc_dir
    type: Directory[]
    label: "Directory of FastQC output files, mirroring input directory structure"
  scanpy_qc_results:
    outputSource: compute_qc_results/scanpy_qc_results
    type: File
    label: "Quality control metrics from Scanpy"
  qc_report:
    outputSource: compute_qc_results/qc_metrics
    type: File
    label: "Quality control report in JSON format"
  dispersion_plot:
    outputSource: scanpy_analysis/dispersion_plot
    type: File
    label: "Gene expression dispersion plot"
  umap_plot:
    outputSource: scanpy_analysis/umap_plot
    type: File
    label: "UMAP dimensionality reduction plot"
  umap_density_plot:
    outputSource: scanpy_analysis/umap_density_plot
    type: File
    label: "UMAP dimensionality reduction plot, colored by cell density"
  spatial_plot:
    outputSource: scanpy_analysis/spatial_plot
    type: File?
    label: "Slide-seq bead plot, colored by Leiden cluster"
  filtered_data_h5ad:
    outputSource: scanpy_analysis/filtered_data_h5ad
    type: File
    label: Full data set of filtered results
    doc: >-
      Full data set of filtered results: expression matrix, coordinates in
      dimensionality-reduced space (PCA and UMAP), cluster assignments via
      the Leiden algorithm, and marker genes for one cluster vs. rest
  marker_gene_plot_t_test:
    outputSource: scanpy_analysis/marker_gene_plot_t_test
    type: File
    label: "Cluster marker genes, t-test"
  marker_gene_plot_logreg:
    outputSource: scanpy_analysis/marker_gene_plot_logreg
    type: File
    label: "Cluster marker genes, logreg method"
  ome_tiff_file
    outputSource: ome_tiff/ome_tiff_file
steps:
  adjust_barcodes:
    in:
      fastq_dir:
        source: fastq_dir
      assay:
        source: assay
    out: [adj_fastq_dir]
    run: steps/adjust-barcodes.cwl
  quantification:
    in:
      fastq_dir:
        source: adjust_barcodes/adj_fastq_dir
      assay:
        source: assay
      threads:
        source: threads
      expected_cell_count:
        source: expected_cell_count
      keep_all_barcodes:
        source: keep_all_barcodes
      visium_probe_set_version:
        source: visium_probe_set_version
    out:
      - salmon_output
      - count_matrix_h5ad
      - raw_count_matrix
      - genome_build_json
    run: steps/salmon-quantification.cwl
  fastqc:
    scatter: [fastq_dir]
    scatterMethod: dotproduct
    in:
      fastq_dir:
        source: fastq_dir
      threads:
        source: threads
    out:
      - fastqc_dir
    run: steps/fastqc.cwl
    label: "Run fastqc on all fastq files in fastq directory"
  scanpy_analysis:
    in:
      assay:
        source: assay
      h5ad_file:
        source: salmon_quantification/count_matrix_h5ad
    out:
      - filtered_data_h5ad
      - umap_plot
      - marker_gene_plot_t_test
      - marker_gene_plot_logreg
      - dispersion_plot
      - umap_density_plot
      - spatial_plot
    run: steps/scanpy-analysis.cwl
    label: "Secondary analysis via ScanPy"
  scvelo_analysis:
    in:
      spliced_h5ad_file:
        source: salmon_quantification/count_matrix_h5ad
      assay_name:
        source: assay
    out:
      - annotated_h5ad_file
      - embedding_grid_plot
    run: steps/scvelo-analysis.cwl
    label: "RNA velocity analysis via scVelo"
  compute_qc_results:
    in:
      assay:
        source: assay
      h5ad_primary:
        source: salmon_quantification/count_matrix_h5ad
      h5ad_secondary:
        source: scanpy_analysis/filtered_data_h5ad
      salmon_dir:
        source: salmon_quantification/salmon_output
    out:
      - scanpy_qc_results
      - qc_metrics
    run: steps/compute-qc-metrics.cwl
    label: "Compute QC metrics"
  ome_tiff:
    in:
      data_dir:
        source: fastq_dir
    out:
      - ome_tiff_file
    run steps/ome-tiff-convert.cwl
    label: "Convert tiff image to ome tiff"
