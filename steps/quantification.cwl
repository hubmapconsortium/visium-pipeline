cwlVersion: v1.0
class: CommandLineTool
requirements:
  DockerRequirement:
    dockerPull: hubmap/visium-quantification
  ResourceRequirement:
    ramMin: 28672
baseCommand: /opt/quantification.py
label: Run Salmon Alevin tool on FASTQ input

# arguments are hardcoded in salmon_wrapper.py

inputs:
  assay:
    type: string
    inputBinding:
      position: 0
  adj_fastq_dir:
    type: Directory
    inputBinding:
      position: 1
  metadata_dir:
    type: Directory
    inputBinding:
      position: 2
  threads:
    type: int
    inputBinding:
      position: 3
      prefix: "--threads"
  expected_cell_count:
    type: int?
    inputBinding:
      position: 4
      prefix: "--expected-cell-count"
  keep_all_barcodes:
    type: boolean?
    inputBinding:
      position: 5
      prefix: "--keep-all-barcodes"
  visium_probe_set_version:
    type: int?
    inputBinding:
      position: 6
      prefix: "--visium-probe-set-version"

outputs:
  h5ad_file:
    type: File?
    outputBinding:
      glob: 'expr.h5ad'
  bam_file:
    type: File?
    outputBinding:
      glob: 'sorted.bam'
  bam_file_index:
    type: File?
    outputBinding:
      glob: 'sorted.bam.bai'
  tsv_gz_file:
    type: File?
    outputBinding:
      glob: 'counts.tsv.gz'
