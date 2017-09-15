#!/bin/bash
# Usage:
#  Run a single conformance test:
#    bash scripts/test_conformance_cwl.bash -n10
#  Run all conformance tests:
#    bash scripts/test_conformance_cwl.bash
#
# This script installs planemo and a Galaxy/planemo cwl-runner implementation
# in the current virtualenv, clones down CWL spec, and runs a configuration test.
# Very few tests will pass - tests which don't have actual data (specify fake jobs)
# will not work, nor will tests which depend on the cwl-runner preserving filenames
# instead of output IDs (this isn't implemented in Galaxy/planemo and I'm not convinced
# it should be), finally any tests of workflows at all or tools using features not
# yet implemented in the CWL will not pass.

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
./run_test.sh $@
