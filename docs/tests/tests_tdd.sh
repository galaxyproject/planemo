#!/bin/bash

set -o xtrace

set -e

# Preconditions - bwa directory doesn't exist and in home directory.
cd
rm -rf bwa

# Tests
conda install --force --yes -c bioconda bwa

echo "Executing bwa to see help - should see output, but it has a non-zero exit code."
bwa | true

planemo project_init --template bwa bwa
cd bwa

cd test-data
bwa index -a is bwa-mem-mt-genome.fa
bwa mem bwa-mem-mt-genome.fa bwa-mem-fastq1.fq bwa-mem-fastq2.fq | \
  samtools view -Sb - > temporary_bam_file.bam && \
  (samtools sort -f temporary_bam_file.bam bwa-aln-test2.bam || samtools sort -o bwa-aln-test2.bam temporary_bam_file.bam)

cd ..

wget https://raw.githubusercontent.com/galaxyproject/planemo/master/docs/writing/bwa-mem_v2.xml -O bwa-mem.xml

echo "Running test we expect to fail..."

planemo t || true

echo "Fixing tool"

wget https://raw.githubusercontent.com/galaxyproject/planemo/master/docs/writing/bwa-mem_v3.xml -O bwa-mem.xml

echo "Running test we expect to be fixed."

planemo t --failed 

# TODO: Finish remaining BWA examples if time.
