"""Unit coverage for the managed PostgreSQL Singularity lifecycle."""

import json
import os
import signal
import subprocess
import sys
from types import SimpleNamespace
from unittest import mock

import click
import pytest
from click.testing import CliRunner

from planemo import options
from planemo.cli import planemo
from planemo.database.factory import started_database_source
from planemo.database.interface import DatabaseConfigurationError
from planemo.database.postgres_singularity import (
    CONTAINER_SOCKET_DIRECTORY,
    DEFAULT_DOCKERIMAGE,
    POSTGRES_SOCKET_NAME,
    SingularityPostgresDatabaseSource,
    start_postgres_singularity,
)
from planemo.io import TERMINATION_TIMEOUT_ENVIRON_KEY
from .test_utils import (
    create_test_context,
    sigterm_ignoring_group,
)


def _source(tmp_path, **kwds):
    source = SingularityPostgresDatabaseSource(
        postgres_storage_location=str(tmp_path / "postgres"),
        **kwds,
    )
    source.startup_timeout = 5
    return source


@pytest.fixture
def singularity_stub(tmp_path):
    """Provide an executable stand-in that records the Singularity CLI contract."""
    calls_path = tmp_path / "singularity-calls.jsonl"
    release_path = tmp_path / "release-container"
    executable = tmp_path / "singularity"
    executable.write_text(f"""#!{sys.executable}
import json
import os
import sys
import time

with open({str(calls_path)!r}, "a") as calls:
    calls.write(json.dumps(sys.argv[1:]) + "\\n")
if "--list" in sys.argv:
    print("postgres | galaxy")
    print("test1234 | galaxy")
elif "run" in sys.argv:
    print("stub container stderr", file=sys.stderr, flush=True)
    print("stub container stdout", flush=True)
    while not os.path.exists({str(release_path)!r}):
        time.sleep(0.01)
""")
    executable.chmod(0o755)

    def calls():
        if not calls_path.exists():
            return []
        return [json.loads(line) for line in calls_path.read_text().splitlines()]

    return SimpleNamespace(command=str(executable), calls=calls, release=release_path)


@click.command()
@options.profile_database_options()
def _database_options_command(**kwds):
    click.echo(f"{kwds['postgres_storage_location']}|{kwds['singularity_cmd']}")


@pytest.mark.parametrize("storage_option", ("--postgres-storage-location", "--postgres_storage_location"))
def test_profile_database_options_accept_storage_aliases_and_singularity_command(tmp_path, storage_option):
    storage = str(tmp_path / "postgres")
    result = CliRunner().invoke(
        _database_options_command,
        [storage_option, storage, "--singularity_cmd", "apptainer"],
        obj=create_test_context(),
    )
    assert result.exit_code == 0, result.output
    assert result.output.strip() == f"{storage}|apptainer"


def test_database_administration_requires_persistent_storage():
    with mock.patch("planemo.database.postgres_singularity.mkdtemp") as make_temp_directory:
        with pytest.raises(DatabaseConfigurationError, match="--postgres-storage-location"):
            started_database_source(database_type="postgres_singularity", for_database_commands=True)
    make_temp_directory.assert_not_called()


def test_database_command_reports_missing_persistent_storage():
    result = CliRunner().invoke(planemo, ["database_list", "--database_type", "postgres_singularity"])
    assert result.exit_code == 2
    assert "requires --postgres-storage-location" in result.output


def test_profile_options_are_owned_by_the_singularity_backend(tmp_path):
    source = _source(tmp_path, singularity_cmd="apptainer", singularity_sudo=False)

    assert source.profile_options() == {
        "postgres_storage_location": str(tmp_path / "postgres"),
        "singularity_cmd": "apptainer",
        "singularity_sudo": False,
    }


def test_container_command_uses_persistent_mounts_environment_and_log(tmp_path, singularity_stub):
    storage = tmp_path / "postgres"
    process = start_postgres_singularity([singularity_stub.command], str(storage))

    try:
        assert os.getpgid(process.pid) == process.pid
    finally:
        singularity_stub.release.touch()
    assert process.wait(timeout=5) == 0
    command = singularity_stub.calls()[0]
    assert command[:2] == ["run", "-B"]
    assert f"{storage / 'pgdata'}:/var/lib/postgresql/data" in command
    assert f"{storage / 'pgrun'}:/var/run/postgresql" in command
    assert "POSTGRES_DB=galaxy" in command
    assert "POSTGRES_USER=galaxy" in command
    assert "POSTGRES_PASSWORD=mysecretpassword" in command
    assert "POSTGRES_INITDB_ARGS=--encoding=UTF-8" in command
    assert f"docker://{DEFAULT_DOCKERIMAGE}" in command
    assert (storage / "postgres.log").read_text().splitlines() == [
        "stub container stderr",
        "stub container stdout",
    ]


def test_start_waits_for_pg_isready_not_just_initialized_cluster(tmp_path):
    source = _source(tmp_path)
    pgdata = tmp_path / "postgres" / "pgdata"
    pgdata.mkdir(parents=True)
    (pgdata / "PG_VERSION").write_text("14")
    process = mock.Mock(pid=42)
    process.poll.return_value = None

    with (
        mock.patch("planemo.database.postgres_singularity.start_postgres_singularity", return_value=process),
        mock.patch.object(source, "_database_is_ready", side_effect=(False, True)) as database_is_ready,
        mock.patch("planemo.database.postgres_singularity.time.sleep") as sleep,
    ):
        source.start()

    assert database_is_ready.call_count == 2
    sleep.assert_called_once_with(1)


def test_readiness_probe_uses_containerized_pg_isready_over_socket(tmp_path, singularity_stub):
    source = _source(tmp_path, singularity_cmd=singularity_stub.command)
    os.makedirs(source.database_socket_dir)
    open(os.path.join(source.database_socket_dir, POSTGRES_SOCKET_NAME), "w").close()

    assert source._database_is_ready()

    command = singularity_stub.calls()[0]
    assert command[:2] == ["exec", "-B"]
    assert f"{source.database_socket_dir}:{CONTAINER_SOCKET_DIRECTORY}" in command
    assert f"docker://{DEFAULT_DOCKERIMAGE}" in command
    assert command[-7:] == [
        "pg_isready",
        "--host",
        CONTAINER_SOCKET_DIRECTORY,
        "--username",
        "galaxy",
        "--dbname",
        "postgres",
    ]


def test_start_reports_container_exit(tmp_path):
    source = _source(tmp_path)
    process = mock.Mock(pid=42)
    process.poll.return_value = 17

    with mock.patch("planemo.database.postgres_singularity.start_postgres_singularity", return_value=process):
        with pytest.raises(RuntimeError, match="code 17"):
            source.start()

    assert source.running_process is None


def test_startup_timeout_stops_the_container(tmp_path):
    source = _source(tmp_path)
    source.startup_timeout = 0
    process = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(300)"], start_new_session=True)
    try:
        with (
            mock.patch("planemo.database.postgres_singularity.start_postgres_singularity", return_value=process),
            mock.patch.object(source, "_database_is_ready", return_value=False),
        ):
            with pytest.raises(RuntimeError, match="did not become ready"):
                source.start()

        assert process.returncode == -signal.SIGTERM
        assert source.running_process is None
    finally:
        if process.poll() is None:
            process.kill()
            process.wait()


def test_stop_waits_then_escalates_the_owned_process_group(tmp_path, monkeypatch):
    source = _source(tmp_path)
    monkeypatch.setenv(TERMINATION_TIMEOUT_ENVIRON_KEY, "0.2")
    with sigterm_ignoring_group(tmp_path / "ready") as process:
        source.running_process = process
        source.stop()

    assert process.returncode == -signal.SIGKILL
    assert source.running_process is None


def test_create_list_and_delete_target_named_databases(tmp_path, singularity_stub):
    source = _source(tmp_path, singularity_cmd=singularity_stub.command)
    storage = tmp_path / "postgres"
    storage.mkdir()
    marker = storage / "must-not-be-deleted"
    marker.write_text("persistent cluster")
    source.create_database("test1234")
    assert source.list_databases() == ["postgres", "test1234"]
    source.delete_database("test1234")

    commands = singularity_stub.calls()
    assert commands[0][-2:] == ["--command", "create database test1234;"]
    assert commands[1][-1] == "--list"
    assert commands[2][-2:] == ["--command", "drop database test1234;"]
    assert marker.read_text() == "persistent cluster"
    for command in commands:
        assert command[0] == "exec"
        assert ["--dbname", "postgres"] == command[command.index("--dbname") : command.index("--dbname") + 2]
