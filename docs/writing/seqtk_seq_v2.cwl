#!/usr/bin/env cwl-runner
cwlVersion: 'v1.0'
class: CommandLineTool
id: "seqtk_seq"
label: "Convert to FASTA (seqtk)"
inputs:
  input1:
    type: File
    doc: |
      TODO
    inputBinding:
      position: 1
      prefix: "-a"
outputs:
  output1:
    type: File
    outputBinding:
      glob: out
baseCommand:
  - "seqtk"
  - "seq"
arguments: []
stdout: out
doc: |
   TODO: Fill in description.