#!/bin/bash
# Usage:
#  Run a single conformance test:
#    bash scripts/test_conformance_cwl.bash -n10
#  Run all conformance tests:
#    bash scripts/test_conformance_cwl.bash
#

set -e
# Ensure working directory is planemo project.
SCRIPTS_DIRECTORY="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_DIRECTORY="${SCRIPTS_DIRECTORY}/.."
cd $PROJECT_DIRECTORY
make install  # will install cwl-runner
pip install cwltest
if [ ! -e common-workflow-language ];
then
    git clone https://github.com/common-workflow-language/common-workflow-language
fi
cd common-workflow-language
git pull
export PLANEMO_ARGS="--no_dependency_resolution"
./run_test.sh --junit-xml=result.xml $@
