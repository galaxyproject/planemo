""" Utilities for using virtualenv as library and planemo command.
"""
from __future__ import absolute_import

import os
import sys

import virtualenv
from galaxy.tools.deps.commands import which


def create_and_exit(virtualenv_path):
    sys.argv = ["virtualenv", virtualenv_path]
    return virtualenv.main()


def create_command(virtualenv_path):
    """ If virtualenv is on Planemo's path use it, otherwise use the planemo
    subcommand virtualenv to create the virtualenv.
    """
    planemo_path = os.path.abspath(sys.argv[0])
    virtualenv_on_path = which("virtualenv")
    if virtualenv_on_path:
        return " ".join([os.path.abspath(virtualenv_on_path), virtualenv_path])
    else:
        return " ".join([planemo_path, "virtualenv", virtualenv_path])
