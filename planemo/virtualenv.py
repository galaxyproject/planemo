"""Utilities for using virtualenv as library and planemo command."""

import os
import sys
from typing import (
    Optional,
    Tuple,
)

from galaxy.util.commands import which

from planemo.io import warn

GALAXY_PYTHON_VERSION_CHOICES = ("3.8", "3.9", "3.10", "3.11", "3.12", "3.13", "3.14")
DEFAULT_PYTHON_VERSION = os.environ.get("PLANEMO_DEFAULT_PYTHON_VERSION", GALAXY_PYTHON_VERSION_CHOICES[-1])


def resolve_python(galaxy_python_version: Optional[str] = None) -> Tuple[str, str, bool]:
    """Return the interpreter used for a Galaxy virtualenv, its version, and
    whether it is a fallback to Planemo's own interpreter.

    If ``pythonX.Y`` is not on the path, Planemo's own interpreter is used.
    """
    if galaxy_python_version is None:
        galaxy_python_version = DEFAULT_PYTHON_VERSION
    python = which("python%s" % galaxy_python_version)
    if python:
        return os.path.abspath(python), galaxy_python_version, False
    return sys.executable or "python", "%d.%d" % sys.version_info[:2], True


def create_command(virtualenv_path: str, galaxy_python_version: Optional[str] = None) -> str:
    """If virtualenv is on Planemo's path use it, otherwise use the planemo
    subcommand virtualenv to create the virtualenv.
    """
    # Create a virtualenv with the selected python version.
    python, python_version, fallback = resolve_python(galaxy_python_version)
    if fallback:
        warn(
            "python%s was not found on PATH, a new Galaxy virtualenv will be created with Planemo's own "
            "interpreter (python %s) instead. Install it or pass --galaxy_python_version to select an "
            "available interpreter.",
            galaxy_python_version or DEFAULT_PYTHON_VERSION,
            python_version,
        )
    virtualenv_on_path = which("virtualenv")
    if virtualenv_on_path:
        command = [virtualenv_on_path, virtualenv_path, "-p", python]
    else:
        command = [python, "-m", "venv", virtualenv_path]
    return " ".join(command)
