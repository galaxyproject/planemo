#!/bin/bash

set -o xtrace

set -e

# Preconditions:
# - planemo installed + toil
# - tool_init_exercise directory doesn't exist and in home directory.
# - conda on PATH

cd
rm -rf tool_init_exercise

# Tests
conda install --force --yes -c bioconda seqtk=1.2

seqtk seq

mkdir tool_init_exercise

cd tool_init_exercise

wget https://raw.githubusercontent.com/galaxyproject/galaxy-test-data/master/2.fastq

seqtk seq -A 2.fastq > 2.fasta

cat 2.fasta

planemo tool_init --cwl --id 'seqtk_seq' --name 'Convert to FASTA (seqtk)'

planemo tool_init --force \
                        --cwl \
                        --id 'seqtk_seq' \
                        --name 'Convert to FASTA (seqtk)' \
                        --example_command 'seqtk seq -A 2.fastq > 2.fasta' \
                        --example_input 2.fastq \
                        --example_output 2.fasta

planemo tool_init --force \
                        --cwl \
                        --id 'seqtk_seq' \
                        --name 'Convert to FASTA (seqtk)' \
                        --example_command 'seqtk seq -A 2.fastq > 2.fasta' \
                        --example_input 2.fastq \
                        --example_output 2.fasta \
                        --requirement seqtk@1.2 \
                        --container 'quay.io/biocontainers/seqtk:1.2--1' \
                        --test_case \
                        --help_from_command 'seqtk seq'

cwltool seqtk_seq.cwl seqtk_seq_job.yml
cwltoil seqtk_seq.cwl seqtk_seq_job.yml

planemo test --no-container seqtk_seq.cwl
planemo test seqtk_seq.cwl
planemo test --no-container --engine toil seqtk_seq.cwl
