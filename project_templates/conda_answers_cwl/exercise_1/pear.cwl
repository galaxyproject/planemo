#!/usr/bin/env cwl-runner
cwlVersion: 'v1.0'
class: CommandLineTool
id: "pear"
label: "Paired-End read merger"

hints:
  SoftwareRequirement:
    packages:
    - package: pear
      version:
      - "0.9.6"

inputs:
  input_f:
    type: File
    doc: |
      Forward reads.
    inputBinding:
      position: 1
      prefix: "-f"
  input_r:
    type: File
    doc: |
      Reverse reads.
    inputBinding:
      position: 2
      prefix: "-r"

  p_value:
    type: float
    default: 0.01
    inputBinding:
      prefix: "--p-value"

  min_overlap:
    type: int
    default: 10
    inputBinding:
      prefix: "--min-overla"

  min_asm_length:
    type: int
    default: 50
    inputBinding:
      prefix: "--min-asm-length"

  min_trim_length:
    type: int
    default: 1
    inputBinding:
      prefix: "--min-trim-length"

  quality_theshold:
    type: int
    default: 0
    inputBinding:
      prefix: "--quality-theshold"

  max_uncalled_base:
    type: float
    default: 1.0
    inputBinding:
      prefix: "--max-uncalled-base"

  test_method:
    type: int
    default: 1
    inputBinding:
      prefix: "--test-method"

  score_method:
    type: int
    default: 2
    inputBinding:
      prefix: "--score-method"

  cap:
    type: int
    default: 40
    inputBinding:
      prefix: "--cap"

outputs:
  assembled_pairs:
    type: File
    outputBinding:
      glob: output.assembled.fastq
  unassembled_forward_reads:
    type: File
    outputBinding:
      glob: output.unassembled.forward.fastq
  unassembled_reverse_reads:
    type: File
    outputBinding:
      glob: output.unassembled.reverse.fastq
  discarded_reads:
    type: File
    outputBinding:
      glob: output.discarded.fastq

baseCommand:
  - "pear"
arguments: ["--phred-base", "33", "-o", "output"]
stdout: out
doc: |
  PEAR is an ultrafast, memory-efficient and highly accurate pair-end read merger.
  PEAR evaluates all possible paired-end read overlaps and without requiring the target fragment
  size as input. In addition, it implements a statistical test for minimizing false-positive results.
  Together with a highly optimized implementation, it can merge millions of paired end reads within
  a couple of minutes on a standard desktop computer.

$namespaces:
  s: http://schema.org/
$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:mainEntity:
  class: s:SoftwareSourceCode
  s:name: PEAR
  s:about: 'Pair-End AssembeR'
  s:url: https://github.com/stevschmid/PEAR
  s:codeRepository: https://github.com/stevschmid/PEAR

  s:license:
  - https://creativecommons.org/licenses/by-sa/3.0/legalcode

  s:targetProduct:
    class: s:SoftwareApplication
    s:softwareVersion: "1.2"
    s:applicationCategory: commandline tool
  s:programmingLanguage: C
 
  s:author:
  - class: s:Person
    s:name: "Jiajie Zhang"
  - class: s:Person
    s:name: "Kassian Kobert"
  - class: s:Person
    s:name: "Tomas Flouri"
  - class: s:Person
    s:name: "Alexandros Stamatakis"

s:downloadUrl: https://github.com/galaxyproject/planemo/blob/master/project_templates/conda_exercises_cwl/pear.cwl
s:codeRepository: https://github.com/galaxyproject/planemo
s:isPartOf:
  class: s:CreativeWork
  s:name: Common Workflow Language
  s:url: http://commonwl.org/

s:author:
  class: s:Person
  s:name: John Chilton
  s:email: mailto:jmchilton@gmail.com
  s:sameAs:
  - id: https://orcid.org/0000-0002-6794-0756
