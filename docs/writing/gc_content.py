#!/usr/bin/env python2
import sys

# usage : python gc_content.py <FASTA file> <output file>

name = None
with open(sys.argv[1]) as infile:
  with open(sys.argv[2], 'w') as outfile:
    for line_raw in infile:
      line = line_raw.rstrip('\r\n')
      if line.startswith('>'):
        if name is not None:
          outfile.write(name+': '+str(gc/float(total))+'\n')
        name = line[1:]
        gc = 0
        total = 0
      else:
        gc += line.count('G') + line.count('C')
        total += len(line)
    outfile.write(name+': '+str(gc/float(total))+'\n')
