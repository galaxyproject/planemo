
``travis_init`` command
======================================

This section is auto-generated from the help text for the planemo command
``travis_init``. This help message can be generated with ``planemo travis_init
--help``.

**Usage**::

    planemo travis_init [OPTIONS] PROJECT

**Help**

Create files to use GitHub/TravisCI testing.

Setup files in a github tool repository to enable continuous
integration testing.

::

    % planemo travis_init .
    % # setup Ubuntu 12.04 w/ dependencies in
    % vim .travis/setup_custom_dependencies.bash
    % git add .travis.yml .travis
    % git commit -m "Add Travis CI testing infrastructure for tools."
    % git push # and register repository @ http://travis-ci.org/

These tests were inspired by work original done and documented by Peter
Cock here http://bit.ly/gxtravisci.

**Options**::


      --help  Show this message and exit.
    
