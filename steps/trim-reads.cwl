cwlVersion: v1.0
class: CommandLineTool
requirements:
  DockerRequirement:
    dockerPull: hubmap/visium-trim-reads:1.0.0
baseCommand: /opt/trim_reads.py
label: Trim FASTQ files

# arguments are hardcoded in salmon_wrapper.py

inputs:
  adj_fastq_dir:
    type: Directory
    inputBinding:
      position: 1
  threads:
    type: int
    inputBinding:
      position: 3
      prefix: "--threads"

outputs:
  trimmed_fastq_dir:
    type: Directory
    outputBinding:
      glob: trimmed
