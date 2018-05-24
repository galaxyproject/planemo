#!/bin/bash

set -o xtrace

set -e

shopt -s expand_aliases  # Needed for conda_env test

# Preconditions:
# - seqtk not on PATH
# - seqtk_example doesn't exist and in home directory
# - conda_exercises doesn't exist in home directory.
# - seqtk environment is absent
# - conda on the PATH

cd

rm -rf seqtk_example
rm -rf conda_exercises
echo "Remove existing seqtk environment - it is fine if it doesn't exist."
conda remove --force --yes --name '__seqtk@1.2' --all || true

echo "Setup completed seqtk example"
# TODO: project_init once in master branch
planemo project_init --template=seqtk_complete_cwl seqtk_example
cd seqtk_example

echo "Check conda_requirements - should see them linting properly"
planemo lint --conda_requirements seqtk_seq.cwl

planemo conda_install seqtk_seq.cwl

echo "Checking for seqtk - shouldn't be present"
which seqtk || true

echo "Simulating '. <(planemo conda_env seqtk_seq.cwl)', subshell-ism doesn't seem to work in script"
planemo conda_env seqtk_seq.cwl > .conda_env
. .conda_env

echo "Inside conda_env I should see seqtk"
which seqtk

conda_env_deactivate

sleep 5  # Need to wait for the package to get cleaned up.

planemo test --no-container seqtk_seq.cwl

cd

# TODO: project_init once in master branch
planemo project_init --template conda_exercises_cwl conda_exercises
cd conda_exercises/exercise_1

ls

echo "Running test wrong - this should fail"
planemo test --no-container pear.cwl || true
echo "Linting pear without requirements - this should fail"
planemo lint --conda_requirements pear.cwl || true

echo "Should see pear from bioconda at 0.9.6"
planemo conda_search pear

wget https://raw.githubusercontent.com/galaxyproject/planemo/master/project_templates_cwl/conda_answers/exercise_1/pear.xml -O pear.xml
planemo lint --conda_requirements pear.cwl || true

planemo test --no-container pear.cwl

cd ../exercise_2

# Fetch and build recipe.

mkdir fleeqtk
cd fleeqtk

wget https://raw.githubusercontent.com/galaxyproject/planemo/master/project_templates/conda_answers_cwl/exercise_2/fleeqtk/build.sh
wget https://raw.githubusercontent.com/galaxyproject/planemo/master/project_templates/conda_answers_cwl/exercise_2/fleeqtk/meta.yaml

cd ..
conda build fleeqtk

# Update to fix tool and see it work.

wget https://raw.githubusercontent.com/galaxyproject/planemo/master/project_templates/conda_answers_cwl/exercise_2/fleeqtk_seq.cwl -O fleeqtk_seq.cwl
planemo conda_install --conda_use_local fleeqtk_seq.cwl
planemo test --no-container fleeqtk_seq.cwl
