cwlVersion: v1.2
class: CommandLineTool
label: Runs NCBI's human read removal tool on each fastq file in fastq directory
requirements:
  DockerRequirement:
    dockerPull: hubmap/human-read-removal:1.1.6
baseCommand: /opt/human_read_removal.py

inputs:
  fastq_dir:
    type: Directory
    doc: Directory containing fastq files to be evaluated
    inputBinding:
      position: 1
  threads:
    type: int
    doc: The number of threads to use for human read removal
    inputBinding:
      position: 2

outputs:
  sanitized_dir:
    type: Directory
    outputBinding:
      glob: "sanitized_fastqs"
    doc: fastq files with any human reads removed
