#!/bin/bash
set -e

planemo=`which planemo`
temp_dir=`mktemp -d`
cd $temp_dir
$planemo conda_init --conda_prefix $temp_dir/conda
export PATH=$temp_dir/conda/bin:$PATH
which conda
conda create -y --name gxtest numpy bx-python pysam
source activate gxtest
git clone --depth 1 https://github.com/galaxyproject/galaxy.git
cd galaxy
./scripts/common_startup.sh --skip-venv --dev-wheels
cd ..
$planemo project_init --template=demo demo
cd demo
$planemo test --galaxy_root ../galaxy --skip_venv cat.xml
