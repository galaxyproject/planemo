#!/usr/bin/env cwl-runner
class: ExpressionTool
requirements:
  - class: InlineJavascriptRequirement
cwlVersion: cwl:draft-3
inputs: []
outputs:
  - { id: output, type: int }
expression: "$({'output': 4})"
