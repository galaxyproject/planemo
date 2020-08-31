""" Utilities for using virtualenv as library and planemo command.
"""
from __future__ import absolute_import

import os
import sys

import virtualenv
from galaxy.util.commands import which


DEFAULT_PYTHON_VERSION = os.environ.get("PLANEMO_DEFAULT_PYTHON_VERSION", "3")


def create_and_exit(virtualenv_path, **kwds):
    sys.argv = ["virtualenv", virtualenv_path]
    python = kwds.get("python", None)
    if python:
        sys.argv.extend(["--python", python])
    return virtualenv.main()


def create_command(virtualenv_path, galaxy_python_version=None):
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

    # Create a virtualenv with the selected python version.
    # default to 2.7
    if galaxy_python_version is None:
        galaxy_python_version = DEFAULT_PYTHON_VERSION
    python = which("python%s" % galaxy_python_version)
    if python:
        python = os.path.abspath(python)
        command.extend(["-p", python])
    command.append(virtualenv_path)
    return " ".join(command)
