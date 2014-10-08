#!/bin/sh
sudo apt-get install -y python-virtualenv
virtualenv planemo-venv
. planemo-venv/bin/activate
pip install git+https://github.com/jmchilton/planemo.git
planemo travis_before_install
 . ${TRAVIS_BUILD_DIR}/.travis/env.sh # source enviornment created by planemo
planemo test --install_galaxy ${TRAVIS_BUILD_DIR}
