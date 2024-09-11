cwlVersion: v1.0
class: CommandLineTool
requirements:
  DockerRequirement:
    dockerPull: hubmap/rna-probes-quantification
  ResourceRequirement:
    ramMin: 28672
baseCommand: /opt/quantification.py
label: Run BWA alignment tool on FASTQ input

# arguments are hardcoded in quantification.py

inputs:
  trimmed_fastq_dir:
    type: Directory
    inputBinding:
      position: 1
  threads:
    type: int
    inputBinding:
      position: 2
      prefix: "--threads"
  organism:
    type: string
    inputBinding:
      position: 3
      prefix: "--organism"

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
