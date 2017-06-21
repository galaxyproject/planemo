#!/bin/bash

set -o xtrace

set -e

# Preconditions - seqtk_example directory doesn't exist and in home directory.
cd
rm -rf seqtk_example
rm -rf conda_testing


echo "Setup completed seqtk example"
planemo project_init --template=seqtk_complete seqtk_example
cd seqtk_example

echo "We should see biocontainer found for this tool"
planemo lint --biocontainers seqtk_seq.xml

echo "This should pass and we should see container was used."
planemo test --biocontainers seqtk_seq.xml

cd ..

planemo project_init --template=conda_testing conda_testing
cd conda_testing/
echo "Should show quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:03dc1d2818d9de56938078b8b78b82d967c1f820-0 is built"
planemo mull bwa_and_samtools.xml


echo "Should see quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:03dc1d2818d9de56938078b8b78b82d967c1f820-0. "
docker images 

echo "Should see cached container was used for this"
planemo test --biocontainers bwa_and_samtools.xml
