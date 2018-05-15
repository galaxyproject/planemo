#!/bin/bash

SCRIPTS_DIRECTORY="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_DIRECTORY="${SCRIPTS_DIRECTORY}/../.."

set -o xtrace

set -e

# Preconditions
# - Planemo installed
# - seqtk_example directory doesn't exist in the home directory.
# - conda_testing directory doesn't exist in the home directory
# - docker running and sudo-less access available via the docker command.

cd
rm -rf conda_exercises
rm -rf conda_testing


echo "Setup completed seqtk example"
# TODO: project_init once in master branch
# planemo project_init --template=conda_exercies_cwl conda_exercies
cp -r $PROJECT_DIRECTORY/project_templates/conda_exercies_cwl conda_exercies
cd conda_exercies

echo "We should see biocontainer found for this tool, but not a Docker container"
planemo lint --biocontainers seqtk_seq.cwl | true

echo "This should pass and we should see container was used."
planemo test --biocontainers --engine cwltool seqtk_seq.cwl

cd ..

planemo project_init --template=conda_testing_cwl conda_testing
cd conda_testing/
echo "Should show quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:03dc1d2818d9de56938078b8b78b82d967c1f820-0 is built"
planemo mull bwa_and_samtools.cwl


echo "Should see quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:03dc1d2818d9de56938078b8b78b82d967c1f820-0. "
docker images 

echo "Should see cached container was used for this"
planemo test --biocontainers bwa_and_samtools.cwl
