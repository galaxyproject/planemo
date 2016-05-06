#!/usr/bin/env cwl-runner
cwlVersion: 'cwl:draft-3'
class: CommandLineTool
id: "seqtk_seq"
label: "Convert to FASTA (seqtk)"
requirements:
  - class: DockerRequirement
    dockerPull: dukegcb/seqtk
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
  