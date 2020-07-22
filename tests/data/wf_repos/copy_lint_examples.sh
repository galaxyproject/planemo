#!/bin/bash

set -x
set -e

PROJECT_DIRECTORY="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
LINT_EXAMPLES="${PROJECT_DIRECTORY}/../../../../gxformat2/tests/examples/lint"

for example_name in "0_basic_format2" "0_basic_native" "1_native_no_output_labels" "1_format2_step_errors";
do
  rm -rf "from_format2/$example_name"
  mkdir "form_format2/$example_name"
  cp "$LINT_EXAMPLES/$example_name.yml" "from_format2/$example_name"
done
