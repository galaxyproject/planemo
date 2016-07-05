#!/bin/sh
sudo apt-get install -y python-virtualenv
virtualenv planemo-venv
. planemo-venv/bin/activate
pip install planemo
planemo travis_before_install
. ${TRAVIS_BUILD_DIR}/.travis/env.sh # source environment created by planemo
planemo test --install_galaxy --no_cache_galaxy ${TRAVIS_BUILD_DIR}
