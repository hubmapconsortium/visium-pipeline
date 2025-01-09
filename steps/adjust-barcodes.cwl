cwlVersion: v1.0
class: CommandLineTool
label: Assay-specific adjustment of cell barcodes
requirements:
  DockerRequirement:
    dockerPull: hubmap/visium-barcode-adj:1.1.2
baseCommand: /opt/adjust_barcodes.py

inputs:
  metadata_dir:
    type: Directory
    inputBinding:
      position: 1
  fastq_dir:
    type: Directory[]
    inputBinding:
      position: 2
outputs:
  adj_fastq_dir:
    type: Directory
    outputBinding:
      glob: 'adj_fastq'
  metadata_json:
    type: File?
    outputBinding:
      glob: 'metadata.json'
