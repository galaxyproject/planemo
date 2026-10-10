"""Unit tests for the Gravity-backed, package-installed Galaxy engine."""

import os
from types import SimpleNamespace
from unittest.mock import (
    Mock,
    patch,
)

import click
import pytest

from planemo.galaxy.config import (
    _stop_daemon_monitor,
    installed_galaxy_config,
    SERVICE_LOG_TAIL_LINES,
    validate_installed_options,
)
from planemo.options import (
    run_engine_option,
    serve_engine_option,
)
from .test_utils import create_test_context


class _Context:
    verbose = False

    def __init__(self, option_sources=None):
        self.option_sources = option_sources or {}

    def get_option_source(self, name, default=None):
        return self.option_sources.get(name, default)

    def vlog(self, message, *args, **kwds):
        pass


@pytest.mark.parametrize(
    ("kwds", "message"),
    [
        ({"galaxy_root": "/checkout"}, "--galaxy_root"),
        ({"install_galaxy": True}, "--install_galaxy"),
        ({"galaxy_branch": "dev"}, "--galaxy_branch"),
        ({"galaxy_branch": "master"}, "--galaxy_branch"),
        ({"galaxy_source": "https://github.com/galaxyproject/galaxy"}, "--galaxy_source"),
    ],
)
def test_installed_option_validation_rejects_checkout_options(kwds, message):
    with pytest.raises(click.UsageError, match=message):
        validate_installed_options(_Context(), kwds)


def test_installed_option_validation_accepts_click_defaults():
    validate_installed_options(
        _Context(),
        {
            "galaxy_root": None,
            "cwl_galaxy_root": None,
            "galaxy_python_version": None,
            "install_galaxy": False,
            "skip_venv": False,
            "no_cache_galaxy": False,
            "galaxy_branch": None,
            "galaxy_source": None,
        },
    )


def test_installed_engine_factory_registration():
    from planemo.engine.factory import (
        build_engine,
        is_galaxy_engine,
    )
    from planemo.engine.galaxy import InstalledGalaxyEngine

    ctx = _Context()

    assert is_galaxy_engine(engine="installed_galaxy")
    assert not is_galaxy_engine(engine="embedded_galaxy")
    assert isinstance(build_engine(ctx, engine="installed_galaxy"), InstalledGalaxyEngine)
    assert isinstance(
        build_engine(ctx, engine="installed_galaxy", database_type="postgres_singularity"), InstalledGalaxyEngine
    )


def test_installed_config_manages_named_database_for_full_lifetime(tmp_path):
    database_source = Mock(keep_running_after_database_commands=False)
    database_source.list_databases.return_value = ["postgres"]
    database_source.sqlalchemy_url.return_value = "postgresql://galaxy@localhost/galaxy"
    config_directory = tmp_path / "config"
    config_directory.mkdir()

    with patch("planemo.database.factory.create_database_source", return_value=database_source):
        with installed_galaxy_config(
            create_test_context(),
            [],
            config_directory=str(config_directory),
            database_type="postgres_singularity",
            port=8765,
        ) as config:
            assert config.galaxy_properties["database_connection"] == "postgresql://galaxy@localhost/galaxy"
            database_source.start.assert_called_once_with()
            database_source.stop.assert_not_called()

    database_source.create_database.assert_called_once_with("galaxy")
    database_source.stop.assert_called_once_with()


@pytest.mark.parametrize("option_factory", [run_engine_option, serve_engine_option])
def test_installed_is_the_only_package_galaxy_engine_choice(option_factory):
    @click.command()
    @option_factory()
    def command(**kwds):
        pass

    engine_option = next(parameter for parameter in command.params if parameter.name == "engine")
    assert "installed_galaxy" in engine_option.type.choices
    assert "embedded_galaxy" not in engine_option.type.choices


def test_missing_gravity_executable_has_actionable_error(tmp_path):
    config_directory = tmp_path / "config"
    config_directory.mkdir()

    with installed_galaxy_config(
        create_test_context(),
        [],
        config_directory=str(config_directory),
        port=8765,
    ) as config:
        with (
            patch("planemo.galaxy.config.sys.executable", str(tmp_path / "bin" / "python")),
            pytest.raises(click.ClickException, match="Gravity's 'galaxy' command"),
        ):
            config.startup_command(_Context())


def test_installed_foreground_serve_uses_managed_process_group(monkeypatch):
    from planemo.galaxy import serve as serve_module

    config = SimpleNamespace(
        env={},
        run_foreground=Mock(return_value=7),
        use_multiprocessing=False,
    )
    monkeypatch.setattr(serve_module, "log_galaxy_command", Mock())

    startup_process, exit_code = serve_module._start_galaxy(_Context(), config, "galaxy command", daemon=False)

    assert startup_process is None
    assert exit_code == 7
    config.run_foreground.assert_called_once_with("galaxy command")


def test_installed_foreground_interruption_cleans_process_group(tmp_path):
    config_directory = tmp_path / "config"
    config_directory.mkdir()
    process = Mock(pid=123)
    process.wait.side_effect = KeyboardInterrupt()

    with installed_galaxy_config(
        create_test_context(),
        [],
        config_directory=str(config_directory),
        port=8765,
    ) as config:
        with (
            patch("planemo.galaxy.config.subprocess.Popen", return_value=process) as popen,
            patch("planemo.galaxy.config.terminate_process_group") as terminate,
            pytest.raises(KeyboardInterrupt),
        ):
            config.run_foreground("galaxy command")

    popen.assert_called_once()
    assert popen.call_args.kwargs["shell"] is True
    assert popen.call_args.kwargs["start_new_session"] is True
    terminate.assert_called_once_with(123, reap=process.poll)


def test_installed_test_config_uses_empty_plugin_directories(tmp_path):
    config_directory = tmp_path / "config"
    config_directory.mkdir()

    with installed_galaxy_config(
        create_test_context(),
        [],
        for_tests=True,
        config_directory=str(config_directory),
        port=8765,
    ) as config:
        properties = config.galaxy_properties
        assert properties["tour_config_dir"] == str(config_directory / "empty")
        assert properties["visualization_plugins_directory"] == str(config_directory / "empty")
        assert properties["interactive_environment_plugins_directory"] == str(config_directory / "empty")


def test_installed_no_cleanup_preserves_config_and_bounded_log_tail():
    lines = [f"installed log line {index}" for index in range(SERVICE_LOG_TAIL_LINES + 5)]

    with installed_galaxy_config(create_test_context(), [], no_cleanup=True, port=8765) as config:
        config_directory = config.config_directory
        with open(config.log_file, "w") as log_fh:
            log_fh.write("\n".join(lines) + "\n")

        service_logs = config.service_log_contents
        assert list(service_logs) == ["installed.log"]
        assert service_logs["installed.log"].splitlines() == lines[-SERVICE_LOG_TAIL_LINES:]

    try:
        assert os.path.isdir(config_directory)
        assert config.log_contents.splitlines() == lines
    finally:
        config.cleanup()


def test_installed_generated_config_is_removed():
    with installed_galaxy_config(create_test_context(), [], port=8765) as config:
        config_directory = config.config_directory
        assert os.path.isdir(config_directory)

    assert not os.path.exists(config_directory)


def test_managed_daemon_monitor_cleanup_is_idempotent(tmp_path, monkeypatch):
    class Process:
        pid = 123
        returncode = None

        def poll(self):
            return self.returncode

    process = Process()
    shutdowns = []

    def shutdown(process_to_stop, asked_to_stop=True):
        shutdowns.append((process_to_stop, asked_to_stop))
        process_to_stop.returncode = 0

    monkeypatch.setattr("planemo.galaxy.config._shut_down_daemon_monitor", shutdown)
    config = SimpleNamespace(
        _daemon_control_fd=None,
        _daemon_process=process,
        pid_file=str(tmp_path / "missing.pid"),
    )

    _stop_daemon_monitor(config)
    _stop_daemon_monitor(config)

    assert shutdowns == [(process, False)]
