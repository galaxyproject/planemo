import os

import click

from planemo.cli import pass_context
from planemo.io import warn, info
from planemo import options
from galaxy.tools.deps.commands import shell

PREPARE_MESSAGE = (
    "Place commands to prepare an Ubuntu VM for use with your tool(s) "
    "in the .travis/setup_custom_dependencies.bash shell script. Be sure to"
    "add these new files to your git repository with 'git add .travis "
    ".travis.yml' and then commit. You will also need to register your github "
    "tool project with Travi CI by visiting https://travis-ci.org/."
)

TRAVIS_YML = """
# This is a special configuration file to run tests on Travis-CI via
# GitHub notifications when changes are committed.
#
# See http://travis-ci.org/ for details

# This file is heavily inspired by Peter Cock's efforts. See more details at
# http://blastedbio.blogspot.com/2013/09/using-travis-ci-for-testing-galaxy-tools.html

# We need Python 2.6 or 2.7 to run Galaxy but many Galaxy tools also need Java
# to run he simplest way to get that installed is to just start with "java".
# (which will still have a system python installed).
language: java

install:
 - sudo apt-get install -y python-virtualenv
 - virtualenv planemo-venv
 - . planemo-venv/bin/activate
 - pip install git+https://github.com/jmchilton/planemo.git
 - planemo travis_before_install
 - . ${TRAVIS_BUILD_DIR}/.travis/env.sh # source enviornment created by planemo

script:
 - planemo test --install_galaxy ${TRAVIS_BUILD_DIR}
"""


@click.command('travis_init')
@options.optional_project_arg()
@pass_context
def cli(ctx, path):
    """Setup files in a github tool repository to enable continuous
    integration testing.::

        % planemo travis_init .
        % # setup Ubuntu 12.04 w/ dependencies in
        % vim .travis/setup_custom_dependencies.bash
        % git add .travis.yml .travis
        % git commit -m "Add Travis CI testing infrastructure for tools."
        % git push # and register repository @ http://travis-ci.org/

    These tests were inspired by work original done and documented by Peter
    Cock here http://bit.ly/gxtravisci.
    """
    shell("mkdir -p '%s/.travis'" % path)
    travis_yml = os.path.join(path, ".travis.yml")
    setup_sh = os.path.join(path, ".travis", "setup_custom_dependencies.bash")
    if not os.path.exists(travis_yml):
        open(travis_yml, "w").write(TRAVIS_YML)
    else:
        warn(".travis.yml file already exists, not overwriting.")
    if not os.path.exists(setup_sh):
        open(setup_sh, "w").write("#!/bin/bash\n")
    info(PREPARE_MESSAGE)
