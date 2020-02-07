"""Module describing the planemo ``travis_init`` command."""
import os

import click

from planemo import options
from planemo import RAW_CONTENT_URL
from planemo.cli import command_function
from planemo.io import info, warn


PREPARE_MESSAGE = (
    "Place commands to prepare an Ubuntu VM for use with your tool(s) "
    "in the .travis/setup_custom_dependencies.bash shell script, which "
    "will be executed on Travis CI via 'planemo travis_before_install'.\n"
    "\n"
    "Be sure to add these new files to your git repository with "
    "'git add .travis .travis.yml' and then commit.\n"
    "\n"
    "You will also need to register your github tool project with "
    "Travis CI by visiting https://travis-ci.org/."
)

TRAVIS_TEST_SCRIPT_URL = RAW_CONTENT_URL + "scripts/travis_test.sh"
TRAVIS_YML = """
# This is a special configuration file to run tests on Travis-CI via
# GitHub notifications when changes are committed.
#
# See http://travis-ci.org/ for details

# This file is heavily inspired by Peter Cock's efforts. See more details at
# http://blastedbio.blogspot.com/2013/09/using-travis-ci-for-testing-galaxy-tools.html

# We need Python 2.6 or 2.7 to run Galaxy but many Galaxy tools also need Java
# to run the simplest way to get that installed is to just start with "java".
# (which will still have a system python installed).
language: java

script:
 - wget -O- %s | /bin/bash -x
""" % TRAVIS_TEST_SCRIPT_URL

TRAVIS_SETUP = """#!/bin/bash
# This will be run on TravisCI via 'planemo travis_before_install' in .travis.yml
#
# TODO: Add your instructions here:
"""


@click.command('travis_init')
@options.optional_project_arg()
@command_function
def cli(ctx, path):
    """Create files to use GitHub/TravisCI testing.

    Setup files in a github tool repository to enable continuous
    integration testing.

    \b
        % planemo travis_init .
        % # setup Ubuntu 12.04 w/ dependencies in
        % vim .travis/setup_custom_dependencies.bash
        % git add .travis.yml .travis
        % git commit -m "Add Travis CI testing infrastructure for tools."
        % git push # and register repository @ http://travis-ci.org/

    These tests were inspired by work original done and documented by Peter
    Cock here http://bit.ly/gxtravisci.
    """
    # TODO: Option --verbose_travis_yaml to unroll travis_test.sh line by line
    # and place all but last in 'install' section and last in 'script'. Would
    # require a yaml dependency though.
    dot_travis_dir = os.path.join(path, '.travis')
    if not os.path.exists(dot_travis_dir):
        os.makedirs(dot_travis_dir)
    travis_yml = os.path.join(path, ".travis.yml")
    setup_sh = os.path.join(dot_travis_dir, "setup_custom_dependencies.bash")
    if not os.path.exists(travis_yml):
        with open(travis_yml, "w") as fh:
            fh.write(TRAVIS_YML)
    else:
        warn(".travis.yml file already exists, not overwriting.")
    if not os.path.exists(setup_sh):
        with open(setup_sh, "w") as fh:
            fh.write(TRAVIS_SETUP)
    else:
        warn("%s already exists, not overwriting." % setup_sh)
    info(PREPARE_MESSAGE)
