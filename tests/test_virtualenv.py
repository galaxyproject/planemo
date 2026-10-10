"""Tests for Planemo's Galaxy virtualenv selection."""

import os
from unittest.mock import (
    call,
    patch,
)

from planemo import virtualenv


def test_galaxy_python_versions_are_explicit_and_current():
    assert virtualenv.GALAXY_PYTHON_VERSION_CHOICES == (
        "3.8",
        "3.9",
        "3.10",
        "3.11",
        "3.12",
        "3.13",
        "3.14",
    )
    assert "3" not in virtualenv.GALAXY_PYTHON_VERSION_CHOICES


def test_default_python_version_uses_latest_supported_version():
    assert virtualenv.DEFAULT_PYTHON_VERSION == os.environ.get("PLANEMO_DEFAULT_PYTHON_VERSION", "3.14")


def test_create_command_uses_default_python_version():
    with (
        patch.object(virtualenv, "DEFAULT_PYTHON_VERSION", "3.14"),
        patch.object(virtualenv, "which", side_effect=["/opt/python3.14", None]) as which,
    ):
        command = virtualenv.create_command("/tmp/galaxy-venv")

    assert command == "/opt/python3.14 -m venv /tmp/galaxy-venv"
    assert which.call_args_list == [call("python3.14"), call("virtualenv")]


def test_create_command_honors_explicit_python_version():
    with patch.object(virtualenv, "which", side_effect=["/opt/python3.10", "/opt/virtualenv"]) as which:
        command = virtualenv.create_command("/tmp/galaxy-venv", "3.10")

    assert command == "/opt/virtualenv /tmp/galaxy-venv -p /opt/python3.10"
    assert which.call_args_list == [call("python3.10"), call("virtualenv")]
