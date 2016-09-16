#!/bin/sh
sudo apt-get install -y python-virtualenv
virtualenv planemo-venv
. planemo-venv/bin/activate
pip install --upgrade pip setuptools
pip install planemo
planemo travis_before_install # runs .travis/setup_custom_dependencies.bash
. ${TRAVIS_BUILD_DIR}/.travis/env.sh # source environment created by planemo
planemo shed_lint --fail_level error --tools --ensure_metadata --urls -r ${TRAVIS_BUILD_DIR}
planemo test --install_galaxy --no_cache_galaxy ${TRAVIS_BUILD_DIR}
