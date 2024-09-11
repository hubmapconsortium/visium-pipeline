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
  metadata_dir:
    label: "Directory containing and metadata.tsv"
    type: Directory
  assay:
    label: "scRNA-seq assay"
    type: string
    default: "10x_v3"
  threads:
    label: "Number of threads for alignment"
    type: int
    default: 1
  organism:
    type: string?
outputs:
  count_matrix_h5ad:
    outputSource: quantification/h5ad_file
    type: File
    label: "Unfiltered count matrix from BWA and umi_tools, converted to H5AD, spliced and unspliced counts"
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
steps:
  adjust_barcodes:
    in:
      metadata_dir:
        source: metadata_dir
      fastq_dir:
        source: fastq_dir
    out: [adj_fastq_dir]
    run: steps/adjust-barcodes.cwl
  trim_reads:
    in:
      adj_fastq_dir:
        source: adjust_barcodes/adj_fastq_dir
      threads:
        source: threads
    out: [trimmed_fastq_dir]
    run: steps/trim-reads.cwl
  quantification:
    in:
      trimmed_fastq_dir:
        source: trim_reads/trimmed_fastq_dir
      threads:
        source: threads
      organism:
        source: organism
    out:
      - h5ad_file
      - bam_file
    run: steps/quantification.cwl
  fastqc:
    scatter: fastq_dir
    scatterMethod: dotproduct
    in:
      fastq_dir:
        source: fastq_dir
      threads:
        source: threads
    out:
      - fastqc_dir
    run: salmon-rnaseq/steps/fastqc.cwl
    label: "Run fastqc on all fastq files in fastq directory"
  scanpy_analysis:
    in:
      assay:
        source: assay
      h5ad_file:
        source: quantification/h5ad_file
    out:
      - filtered_data_h5ad
      - umap_plot
      - marker_gene_plot_t_test
      - marker_gene_plot_logreg
      - dispersion_plot
      - umap_density_plot
      - spatial_plot
    run: salmon-rnaseq/steps/scanpy-analysis.cwl
    label: "Secondary analysis via ScanPy"
  compute_qc_results:
    in:
      h5ad_primary:
        source: quantification/h5ad_file
      h5ad_secondary:
        source: scanpy_analysis/filtered_data_h5ad
      bam_file:
        source: quantification/bam_file
    out:
      - scanpy_qc_results
      - qc_metrics
    run: steps/compute-qc-metrics.cwl
    label: "Compute QC metrics"

