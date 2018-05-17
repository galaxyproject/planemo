#!/usr/bin/env cwl-runner
class: ExpressionTool
requirements:
  - class: InlineJavascriptRequirement
cwlVersion: 'v1.0'
inputs: []
outputs:
  - { id: output, type: int }
expression: "$({'output': 4})"
