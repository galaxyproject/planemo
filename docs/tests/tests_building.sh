#!/bin/bash

set -o xtrace

set -e

# Preconditions:
# - planemo installed
# - tool_init_exercise directory doesn't exist and in home directory.
# - local tool shed running and configured with ~/.planemo.yml
# - conda on PATH

cd
rm -rf tool_init_exercise

# Tests

conda install --force --yes -c conda-forge -c bioconda seqtk=1.2

cat <(seqtk seq)

mkdir tool_init_exercise

cd tool_init_exercise

wget https://raw.githubusercontent.com/galaxyproject/galaxy-test-data/master/2.fastq

seqtk seq -A 2.fastq > 2.fasta

cat 2.fasta

planemo tool_init --id 'seqtk_seq' --name 'Convert to FASTA (seqtk)'

cat seqtk_seq.xml

planemo tool_init --force \
                    --id 'seqtk_seq' \
                    --name 'Convert to FASTA (seqtk)' \
                    --requirement seqtk@1.2 \
                    --example_command 'seqtk seq -a 2.fastq > 2.fasta' \
                    --example_input 2.fastq \
                    --example_output 2.fasta

cat seqtk_seq.xml

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

planemo l

planemo t

planemo shed_init --name=seqtk_seq \
                    --owner=planemo \
                    --description=seqtk_seq \
                    --long_description="Tool that converts FASTQ to FASTA files using seqtk" \
                    --category="Fastq Manipulation"

planemo shed_lint --tools

planemo shed_create --shed_target local

planemo shed_update --shed_target local

planemo tool_init --force \
                    --macros \
                    --id 'seqtk_seq' \
                    --name 'Convert to FASTA (seqtk)' \
                    --requirement seqtk@1.2 \
                    --example_command 'seqtk seq -A 2.fastq > 2.fasta' \
                    --example_input 2.fastq \
                    --example_output 2.fasta \
                    --test_case \
                    --help_from_command 'seqtk seq'

cat seqtk_seq.xml

