#!/bin/bash

set -o xtrace

set -e

# Preconditions
# - seqtk not on PATH
# - Planemo installed
# - conda_exercises directory doesn't exist in the home directory. 
# - conda_answers directory doesn't exist in the home directory.
# - conda_testing directory doesn't exist in the home directory
# - docker running and sudo-less access available via the docker command.

cd
rm -rf conda_exercises
rm -rf conda_answers_cwl
rm -rf conda_testing

echo "Setup completed seqtk example"
planemo project_init --template=conda_exercies_cwl conda_exercies
cd conda_exercies

cd conda_exercises/exercise_3

echo "We should see biocontainer found for this tool, but not a Docker container"
planemo lint --biocontainers seqtk_seq.cwl | true

echo "This should pass and we should see container was used."
planemo test --biocontainers seqtk_seq.cwl

echo "This should fail without biocontainers"
planemo test seqtk_seq.cwl

cd ..

echo "Run seqtk exercise answers"
planemo project_init --template=conda_answers_cwl conda_answers
cd conda_answers/exercise_3

planemo lint --biocontainers seqtk_seq.cwl
planemo test seqtk_seq.cwl


planemo project_init --template=conda_testing_cwl conda_testing
cd conda_testing/
echo "Should show quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:03dc1d2818d9de56938078b8b78b82d967c1f820-0 is built"
planemo mull bwa_and_samtools.cwl


echo "Should see quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:03dc1d2818d9de56938078b8b78b82d967c1f820-0. "
docker images 

echo "Should see cached container was used for this"
planemo test --biocontainers bwa_and_samtools.cwl
