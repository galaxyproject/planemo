"""CLI regressions for selecting the installed Galaxy runtime."""

import importlib
import json
import shutil
from unittest.mock import Mock

import pytest
from click.testing import CliRunner

from planemo.cli import planemo
from .test_utils import PROJECT_TEMPLATES_DIR

autoupdate_command = importlib.import_module("planemo.commands.cmd_autoupdate")
test_command = importlib.import_module("planemo.commands.cmd_test")


@pytest.fixture
def workspace(tmp_path, monkeypatch):
    monkeypatch.setenv("PLANEMO_GLOBAL_CONFIG_PATH", str(tmp_path / "planemo.yml"))
    return tmp_path / "workspace"


@pytest.mark.parametrize("engine", ["galaxy", "installed_galaxy"])
@pytest.mark.parametrize("database_type", ["sqlite", "postgres", "postgres_singularity"])
def test_profile_retains_local_engine(tmp_path, monkeypatch, workspace, engine, database_type):
    source = Mock(store_connection_in_profile=database_type == "postgres")
    source.sqlalchemy_url.return_value = "postgresql://galaxy@localhost/galaxy"
    source.profile_options.return_value = {}
    monkeypatch.setattr("planemo.database.factory.create_database_source", Mock(return_value=source))
    runner = CliRunner()
    result = runner.invoke(
        planemo,
        [
            "--directory",
            str(workspace),
            "profile_create",
            "local",
            "--engine",
            engine,
            "--database_type",
            database_type,
        ],
    )
    assert result.exit_code == 0, result.output
    profile = workspace / "profiles" / "local" / "planemo_profile_options.json"
    assert json.loads(profile.read_text())["engine"] == engine

    run_tests = Mock(return_value=0)
    monkeypatch.setattr(test_command, "test_runnables", run_tests)
    tool = str(PROJECT_TEMPLATES_DIR + "/demo/cat.xml")
    result = runner.invoke(planemo, ["--directory", str(workspace), "test", "--profile", "local", tool])
    assert result.exit_code == 0, result.output
    assert run_tests.call_args.kwargs["engine"] == engine


@pytest.mark.parametrize("engine", [None, "galaxy", "installed_galaxy", "docker_galaxy", "external_galaxy"])
def test_autoupdate_verifies_with_selected_engine(tmp_path, monkeypatch, workspace, engine):
    tool = tmp_path / "cat.xml"
    shutil.copyfile(PROJECT_TEMPLATES_DIR + "/demo/cat.xml", tool)
    monkeypatch.setattr(autoupdate_command.autoupdate, "autoupdate_tool", Mock(return_value=[str(tool)]))
    run_tests = Mock(return_value=0)
    monkeypatch.setattr(autoupdate_command, "test_runnables", run_tests)
    arguments = ["--directory", str(workspace), "autoupdate", "--test", str(tool)]
    if engine:
        arguments.extend(["--engine", engine])
    result = CliRunner().invoke(planemo, arguments)
    assert result.exit_code == 0, result.output
    assert run_tests.call_args.kwargs["engine"] == (engine or "galaxy")
