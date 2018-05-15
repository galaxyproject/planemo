#!/usr/bin/env cwl-runner
cwlVersion: 'v1.0'
class: CommandLineTool
id: "bwa_and_samtools"

doc:
  Dump bwa and samtools default help text to text datasets.

hints:
  SoftwareRequirement:
    packages:
    - package: bwa
      version:
      - "0.7.15"
    - package: samtools
      version:
      - "1.3.1"

inputs: {}
outputs:
outputs:
  bwa_help:
    type: File
    outputBinding:
      glob: bwa_help.txt
  
  samtools_help:
    type: File
    outputBinding:
      glob: samtools_help.txt

successCodes: [0, 1]  # samtool help message returns exit code of 1

baseCommand:
 - sh
 - "-c"
 - "bwa > bwa_help.txt 2>&1; samtools > samtools_help.txt 2>&1"
