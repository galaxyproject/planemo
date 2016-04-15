#!/bin/bash
set -e
# Ensure working directory is planemo project.
SCRIPTS_DIRECTORY="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_DIRECTORY="${SCRIPTS_DIRECTORY}/.."
cd $PROJECT_DIRECTORY
make install  # will install cwl-runner
if [ ! -e common-workflow-language ];
then
    git clone https://github.com/common-workflow-language/common-workflow-language
fi
cd common-workflow-language
git pull
./run_test.sh $@
