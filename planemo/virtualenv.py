""" Utilities for using virtualenv as library and planemo command.
"""
from __future__ import absolute_import

import os
import sys

import virtualenv
from galaxy.tools.deps.commands import which


def create_and_exit(virtualenv_path, **kwds):
    sys.argv = ["virtualenv", virtualenv_path]
    python = kwds.get("python", None)
    if python:
        sys.argv.extend(["--python", python])
    return virtualenv.main()


def create_command(virtualenv_path):
    """ If virtualenv is on Planemo's path use it, otherwise use the planemo
    subcommand virtualenv to create the virtualenv.
    """
    planemo_path = os.path.abspath(sys.argv[0])
    virtualenv_on_path = which("virtualenv")
    if virtualenv_on_path:
        base_command = [
            os.path.abspath(virtualenv_on_path),
        ]
    else:
        base_command = [
            planemo_path, "virtualenv",
        ]

    command = base_command

    # If planemo is running in a Python 3 environment but Python 2.7
    # is available for Galaxy, use it.
    python27 = which("python2.7")
    if python27:
        python27 = os.path.abspath(python27)
        command.extend(["-p", python27])
    command.append(virtualenv_path)
    return " ".join(command)
