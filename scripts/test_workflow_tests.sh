#!/bin/bash

# Ensure working directory is planemo project.
SCRIPTS_DIRECTORY="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_DIRECTORY="${SCRIPTS_DIRECTORY}/.."

TEMP_DIRECTORY=`mktemp -d 2>/dev/null || mktemp -d -t 'planemowftest'`
cp "$SCRIPTS_DIRECTORY/run_galaxy_workflow_tests.sh" "$TEMP_DIRECTORY"
cp "$PROJECT_DIRECTORY/tests/data/"wf4* "$TEMP_DIRECTORY"
cp "$PROJECT_DIRECTORY/tests/data/1.bed" "$TEMP_DIRECTORY"

cd $PROJECT_DIRECTORY

# Build Planemo wheel.
make dist

# Test against wheel.
export PLANEMO_TARGET="$PROJECT_DIRECTORY/dist/planemo*whl"

cd $TEMP_DIRECTORY
bash run_workflow_tests.sh "wf4-distro-tools.gxwf.yml"
