#!/usr/bin/env cwl-runner
cwlVersion: 'cwl:draft-3'
class: CommandLineTool
id: "seqtk_seq"
label: "Convert to FASTA (seqtk)"
inputs:
  - id: input1
    type: File
    description: |
      TODO
    inputBinding:
      position: 1
      prefix: "-a"
outputs:
  - id: output1
    type: File
    outputBinding:
      glob: out
baseCommand:
  - "seqtk"
  - "seq"
arguments: []
stdout: out
description: |
   TODO: Fill in description.