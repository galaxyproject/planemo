#!/bin/bash

# Usage: DOC_TEST_RUNNER="https://raw.githubusercontent.com/galaxyproject/planemo/master/scripts/run_doc_test.sh"
#        DOCS=building bash <(curl -s "$DOC_TEST_RUNNER")
#        DOCS=conda bash <(curl -s "$DOC_TEST_RUNNER")
#        DOCS=conda_cwl bash <(curl -s "$DOC_TEST_RUNNER")

set -e

: ${PLANEMO_TARGET:="."}
: ${PLANEMO_VIRTUAL_ENV:=".venv-doc-tests"}
: ${DOCS:="building"}

# Ensure Planemo is installed.
if [ ! -d "${PLANEMO_VIRTUAL_ENV}" ]; then
    virtualenv "${PLANEMO_VIRTUAL_ENV}"
    . "${PLANEMO_VIRTUAL_ENV}"/bin/activate
    pip install -U pip>7
    # Intentionally expand wildcards in PLANEMO_TARGET.
    shopt -s extglob
    pip install ${PLANEMO_TARGET}
fi
. "${PLANEMO_VIRTUAL_ENV}"/bin/activate

planemo conda_init | true
export PATH="$HOME/miniconda3/bin:$PATH"

bash docs/tests/tests_"${DOCS}".sh
