#!/bin/bash
set -e

temp_dir=`mktemp -d`
cd $temp_dir
planemo conda_init --conda_prefix $temp_dir/conda
export PATH=$temp_dir/conda/bin:$PATH
which conda
conda create --name gxtest
source activate gxtest
conda install -y numpy bx-python pysam
git clone https://github.com/galaxyproject/galaxy.git
cd galaxy
./scripts/common_startup.sh --skip-venv --dev-wheels
cd ..
planemo project_init --template=demo demo
cd demo
planemo test --galaxy_root ../galaxy --skip_venv cat.xml
