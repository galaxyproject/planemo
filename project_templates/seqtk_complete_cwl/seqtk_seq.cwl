#!/usr/bin/env cwl-runner
cwlVersion: 'v1.0'
class: CommandLineTool
id: "seqtk_seq"
label: "Convert to FASTA (seqtk)"
hints:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/seqtk:1.2--1
  SoftwareRequirement:
    packages:
    - package: seqtk
      version:
      - "1.2"
inputs:
  input1:
    type: File
    doc: |
      Input FASTA file.
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
  
  Usage:   seqtk seq [options] <in.fq>|<in.fa>

  Options: -q INT    mask bases with quality lower than INT [0]
           -X INT    mask bases with quality higher than INT [255]
           -n CHAR   masked bases converted to CHAR; 0 for lowercase [0]
           -l INT    number of residues per line; 0 for 2^32-1 [0]
           -Q INT    quality shift: ASCII-INT gives base quality [33]
           -s INT    random seed (effective with -f) [11]
           -f FLOAT  sample FLOAT fraction of sequences [1]
           -M FILE   mask regions in BED or name list FILE [null]
           -L INT    drop sequences with length shorter than INT [0]
           -c        mask complement region (effective with -M)
           -r        reverse complement
           -A        force FASTA output (discard quality)
           -C        drop comments at the header lines
           -N        drop sequences containing ambiguous bases
           -1        output the 2n-1 reads only
           -2        output the 2n reads only
           -V        shift quality by '(-Q) - 33'
           -U        convert all bases to uppercases
           -S        strip of white spaces in sequences

$namespaces:
  s: http://schema.org/
$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:mainEntity:
  class: s:SoftwareSourceCode
  s:name: seqtk
  s:about: 'Toolkit for processing sequences in FASTA/Q formats'
  s:url: https://github.com/lh3/seqtk
  s:codeRepository: https://github.com/lh3/seqtk

  s:license:
  - https://opensource.org/licenses/MIT

  s:targetProduct:
    class: s:SoftwareApplication
    s:softwareVersion: "1.2"
    s:applicationCategory: commandline tool
  s:programmingLanguage: C
  s:publication:
  - class: s:ScholarlyArticle
    id: http://dx.doi.org/10.1093/bioinformatics/bts635

  s:author:
  - class: s:Person
    s:name: "Heng Li"
    s:sameAs:
    - id: https://orcid.org/0000-0003-4874-2874

s:downloadUrl: https://github.com/galaxyproject/planemo/blob/master/project_templates/seqtk_complete_cwl/seqtk_seq.cwl
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
