#!/bin/bash

planemo tool_init --id 'seqtk_seq' --name 'Convert to FASTA (seqtk)'
mv seqtk_seq.xml seqtk_seq_v1.xml
planemo tool_init --force \
                  --id 'seqtk_seq' \
                  --name 'Convert to FASTA (seqtk)' \
                  --requirement seqtk@1.2 \
                  --example_command 'seqtk seq -a 2.fastq > 2.fasta' \
                  --example_input 2.fastq \
                  --example_output 2.fasta
mv seqtk_seq.xml seqtk_seq_v2.xml
planemo tool_init --force \
                  --id 'seqtk_seq' \
                  --name 'Convert to FASTA (seqtk)' \
                  --requirement seqtk@1.2 \
                  --example_command 'seqtk seq -a 2.fastq > 2.fasta' \
                  --example_input 2.fastq \
                  --example_output 2.fasta \
                  --test_case \
                  --cite_url 'https://github.com/lh3/seqtk' \
                  --help_from_command 'seqtk seq'
mv seqtk_seq.xml seqtk_seq_v3.xml


planemo tool_init --cwl --id 'seqtk_seq' --name 'Convert to FASTA (seqtk)'
mv seqtk_seq.cwl seqtk_seq_v1.cwl

planemo tool_init --force \
                  --cwl \
                  --id 'seqtk_seq' \
                  --name 'Convert to FASTA (seqtk)' \
                  --example_command 'seqtk seq -a 2.fastq > 2.fasta' \
                  --example_input 2.fastq \
                  --example_output 2.fasta
mv seqtk_seq.cwl seqtk_seq_v2.cwl

planemo tool_init --force \
                  --cwl \
                  --id 'seqtk_seq' \
                  --name 'Convert to FASTA (seqtk)' \
                  --example_command 'seqtk seq -a 2.fastq > 2.fasta' \
                  --example_input 2.fastq \
                  --example_output 2.fasta \
                  --test_case \
                  --container 'dukegcb/seqtk' \
                  --help_from_command 'seqtk seq'                  
mv seqtk_seq.cwl seqtk_seq_v3.cwl
mv seqtk_seq_tests.yml seqtk_seq_tests_v3.yml
