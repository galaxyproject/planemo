#!/bin/bash

set -o xtrace

set -e

DOC_TESTS_DIRECTORY="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

bash "${DOC_TESTS_DIRECTORY}/tests_building.sh"

bash "${DOC_TESTS_DIRECTORY}/tests_building_cwl.sh"

bash "${DOC_TESTS_DIRECTORY}/tests_tdd.sh"

bash "${DOC_TESTS_DIRECTORY}/tests_conda.sh"

bash "${DOC_TESTS_DIRECTORY}/tests_docker.sh"
