import contextlib
import importlib
import json
import os
import signal
import subprocess
import sys
import time

import pytest
from packaging.version import InvalidVersion

from planemo import network_util
from planemo.galaxy.config import (
    gravity_supports_multiprocessing,
    LocalGalaxyConfig,
    write_galaxy_config,
)
from .test_utils import create_test_context

serve_module = importlib.import_module("planemo.galaxy.serve")


def _make_galaxy_root(tmp_path, galaxy_version="25.0.1"):
    galaxy_root = tmp_path / "galaxy"
    dependencies = galaxy_root / "lib" / "galaxy" / "dependencies"
    dependencies.mkdir(parents=True)
    if galaxy_version is not None:
        version_package = galaxy_root / "lib" / "galaxy" / "version"
        version_package.mkdir()
        version_parts = galaxy_version.split(".", 2)
        version_major = ".".join(version_parts[:2])
        version_minor = version_parts[2] if len(version_parts) == 3 else ""
        (version_package / "__init__.py").write_text(
            f'VERSION_MAJOR = "{version_major}"\n'
            f'VERSION_MINOR = "{version_minor}"\n'
            'VERSION = VERSION_MAJOR + (f".{VERSION_MINOR}" if VERSION_MINOR else "")\n'
        )
    return galaxy_root


@pytest.mark.parametrize(
    ("galaxy_version", "expected"),
    [
        ("24.2.4", False),
        ("25.0.0", False),
        ("25.0.1", True),
        ("25.1.0", True),
    ],
)
def test_gravity_multiprocessing_detection(tmp_path, galaxy_version, expected):
    galaxy_root = _make_galaxy_root(tmp_path, galaxy_version)
    assert gravity_supports_multiprocessing(str(galaxy_root)) is expected


@pytest.mark.parametrize(
    ("galaxy_version", "exception"),
    [("not-a-version", InvalidVersion), (None, FileNotFoundError)],
)
def test_gravity_multiprocessing_detection_fails_for_invalid_galaxy_version(tmp_path, galaxy_version, exception):
    galaxy_root = _make_galaxy_root(tmp_path, galaxy_version)
    with pytest.raises(exception):
        gravity_supports_multiprocessing(str(galaxy_root))


@pytest.mark.parametrize("galaxy_version", ["25.0.0", "25.0.1"])
def test_gravity_configuration_selects_process_manager(tmp_path, galaxy_version):
    galaxy_root = _make_galaxy_root(tmp_path, galaxy_version)
    config_directory = tmp_path / "config"
    config_directory.mkdir()
    env = {"SUPERVISORD_SOCKET": "old-socket"}

    write_galaxy_config(
        str(galaxy_root),
        {},
        env,
        {},
        {"port": 12345},
        lambda name: str(config_directory / name),
    )

    with open(env["GALAXY_CONFIG_FILE"]) as config_file:
        gravity = json.load(config_file)["gravity"]
    if galaxy_version == "25.0.1":
        assert gravity["process_manager"] == "multiprocessing"
        assert "SUPERVISORD_SOCKET" not in env
    else:
        assert "process_manager" not in gravity
        assert env["SUPERVISORD_SOCKET"] != "old-socket"


def _make_local_config(galaxy_root, port):
    config_directory = galaxy_root / "planemo-config"
    config_directory.mkdir()
    env = {
        "GALAXY_LOG": str(galaxy_root / "main.log"),
        "GALAXY_PID": str(galaxy_root / "main.pid"),
        "GRAVITY_STATE_DIR": str(config_directory / "gravity"),
        "TEST_RECORD": str(galaxy_root / "record.json"),
    }
    os.makedirs(env["GRAVITY_STATE_DIR"])
    return LocalGalaxyConfig(
        create_test_context(),
        str(config_directory),
        env,
        None,
        port,
        "main",
        "test-key",
        [],
        str(galaxy_root),
        {},
    )


def _write_executable(path, contents):
    path.write_text(contents)
    path.chmod(0o755)


MODERN_RUN_SH = r"""#!/usr/bin/env python3
import json
import os
import signal
import subprocess
import sys
from http.server import BaseHTTPRequestHandler, HTTPServer

record_path = os.environ["TEST_RECORD"]
child = subprocess.Popen(
    [sys.executable, "-c", 'import time; print("modern child output", flush=True); time.sleep(300)']
)
with open(record_path, "w") as record:
    json.dump({
        "args": sys.argv[1:],
        "child_pid": child.pid,
        "env": {key: os.environ.get(key) for key in ["GALAXY_LOG", "GALAXY_PID", "SUPERVISORD_SOCKET"]},
        "pid": os.getpid(),
    }, record)
print("modern daemon started", flush=True)

def stop(signum, frame):
    try:
        child.wait(timeout=2)
    except subprocess.TimeoutExpired:
        child.kill()
        child.wait()
    raise SystemExit(0)

signal.signal(signal.SIGTERM, stop)

class Handler(BaseHTTPRequestHandler):
    def do_GET(self):
        if self.path == "/api/version":
            body = b'{"version_major": "25.0"}'
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)
        else:
            self.send_error(404)
    def log_message(self, format, *args):
        pass

server = HTTPServer(("127.0.0.1", int(os.environ["TEST_PORT"])), Handler)
server.serve_forever()
"""


STUBBORN_RUN_SH = r"""#!/usr/bin/env python3
import json
import os
import signal
import subprocess
import sys
import time

child = subprocess.Popen(
    [
        sys.executable,
        "-c",
        "import signal, time; signal.signal(signal.SIGTERM, signal.SIG_IGN); time.sleep(300)",
    ]
)
with open(os.environ["TEST_RECORD"], "w") as record:
    json.dump({"child_pid": child.pid, "pid": os.getpid()}, record)
signal.signal(signal.SIGTERM, signal.SIG_IGN)
while True:
    time.sleep(60)
"""


LEGACY_RUN_SH = r"""#!/usr/bin/env python3
import json
import os
import subprocess
import sys

with open(os.environ["TEST_RECORD"], "w") as record:
    json.dump({"args": sys.argv[1:], "env": {"SUPERVISORD_SOCKET": os.environ.get("SUPERVISORD_SOCKET")}}, record)
process = subprocess.Popen(
    [sys.executable, os.path.join(os.path.dirname(__file__), "server.py")],
    env=os.environ.copy(),
    start_new_session=True,
)
with open(os.environ["GALAXY_PID"], "w") as pid_file:
    pid_file.write(str(process.pid))
with open(os.path.join(os.environ["GRAVITY_STATE_DIR"], "legacy.pid"), "w") as pid_file:
    pid_file.write(str(process.pid))
"""


SERVER_SCRIPT = r"""#!/usr/bin/env python3
import os
from http.server import BaseHTTPRequestHandler, HTTPServer

class Handler(BaseHTTPRequestHandler):
    def do_GET(self):
        body = b'{"version_major": "25.0"}'
        self.send_response(200)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)
    def log_message(self, format, *args):
        pass

HTTPServer(("127.0.0.1", int(os.environ["TEST_PORT"])), Handler).serve_forever()
"""


GALAXYCTL = r"""#!/usr/bin/env python3
import json
import os
import signal
import sys

with open(os.environ["TEST_GALAXYCTL_RECORD"], "w") as record:
    json.dump({"args": sys.argv[1:]}, record)
with open(os.path.join(os.environ["GRAVITY_STATE_DIR"], "legacy.pid")) as pid_file:
    os.killpg(int(pid_file.read()), signal.SIGTERM)
"""


@contextlib.contextmanager
def _configured_galaxy(config):
    yield config


def _pid_exists(pid):
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    return True


def _wait_for_pid_exit(pid):
    for _ in range(100):
        if not _pid_exists(pid):
            return
        time.sleep(0.05)
    raise AssertionError(f"Process {pid} did not exit")


def _wait_for_path(path):
    for _ in range(100):
        if path.exists():
            return
        time.sleep(0.05)
    raise AssertionError(f"Path {path} was not created")


def _wait_for_json(path):
    for _ in range(100):
        try:
            return json.loads(path.read_text())
        except (FileNotFoundError, json.JSONDecodeError):
            time.sleep(0.05)
    raise AssertionError(f"Valid JSON was not written to {path}")


def test_multiprocessing_daemon_uses_foreground_process_and_cleans_group(tmp_path, monkeypatch):
    galaxy_root = _make_galaxy_root(tmp_path, "25.0.1")
    _write_executable(galaxy_root / "run.sh", MODERN_RUN_SH)
    port = network_util.get_free_port()
    config = _make_local_config(galaxy_root, port)
    config.env["TEST_PORT"] = str(port)
    external_pid = tmp_path / "external.pid"
    monkeypatch.setattr(serve_module, "galaxy_config", lambda *args, **kwds: _configured_galaxy(config))
    monkeypatch.setattr(importlib.import_module("planemo.galaxy.ephemeris_sleep"), "DEFAULT_SLEEP_WAIT", 0.05)

    with serve_module.serve_daemon(
        config._ctx,
        daemon=True,
        port=port,
        pid_file=str(external_pid),
        skip_venv=True,
        galaxy_startup_timeout=200,
    ):
        with open(config.env["TEST_RECORD"]) as record_file:
            record = json.load(record_file)
        assert "--daemon" not in record["args"]
        assert record["args"] == ["--server-name", "main"]
        assert record["env"]["GALAXY_LOG"] == config.log_file
        assert record["env"]["SUPERVISORD_SOCKET"] is None
        assert "modern daemon started" in config.log_contents
        assert "modern child output" in config.log_contents
        assert config.service_log_contents == {}
        assert os.path.islink(external_pid)
        assert int(external_pid.read_text()) == config._daemon_process.pid
        assert _pid_exists(record["pid"])
        assert _pid_exists(record["child_pid"])

    _wait_for_pid_exit(record["pid"])
    _wait_for_pid_exit(record["child_pid"])
    assert config._daemon_process.returncode is not None
    assert not os.path.exists(config.pid_file)


def test_legacy_daemon_uses_supervisor_lifecycle(tmp_path, monkeypatch):
    galaxy_root = _make_galaxy_root(tmp_path, "25.0.0")
    _write_executable(galaxy_root / "run.sh", LEGACY_RUN_SH)
    _write_executable(galaxy_root / "server.py", SERVER_SCRIPT)
    galaxyctl = galaxy_root / ".venv" / "bin" / "galaxyctl"
    galaxyctl.parent.mkdir(parents=True)
    _write_executable(galaxyctl, GALAXYCTL)
    port = network_util.get_free_port()
    config = _make_local_config(galaxy_root, port)
    config.env.update(
        {
            "SUPERVISORD_SOCKET": str(tmp_path / "supervisor.sock"),
            "TEST_GALAXYCTL_RECORD": str(galaxy_root / "galaxyctl.json"),
            "TEST_PORT": str(port),
        }
    )
    monkeypatch.setattr(serve_module, "galaxy_config", lambda *args, **kwds: _configured_galaxy(config))
    monkeypatch.setattr(importlib.import_module("planemo.galaxy.ephemeris_sleep"), "DEFAULT_SLEEP_WAIT", 0.05)

    with serve_module.serve_daemon(
        config._ctx,
        daemon=True,
        port=port,
        skip_venv=True,
        galaxy_startup_timeout=200,
    ):
        with open(config.env["TEST_RECORD"]) as record_file:
            record = json.load(record_file)
        assert record["args"] == ["--daemon"]
        assert record["env"]["SUPERVISORD_SOCKET"] == config.env["SUPERVISORD_SOCKET"]

    with open(config.env["TEST_GALAXYCTL_RECORD"]) as record_file:
        galaxyctl_record = json.load(record_file)
    assert galaxyctl_record["args"] == ["--state-dir", config.gravity_state_dir, "shutdown"]


def test_multiprocessing_startup_exit_reports_log_without_waiting_for_timeout(tmp_path, monkeypatch):
    galaxy_root = _make_galaxy_root(tmp_path, "25.0.1")
    _write_executable(galaxy_root / "run.sh", "#!/bin/sh\necho controlled startup failure\nexit 7\n")
    port = network_util.get_free_port()
    config = _make_local_config(galaxy_root, port)
    monkeypatch.setattr(serve_module, "galaxy_config", lambda *args, **kwds: _configured_galaxy(config))
    monkeypatch.setattr(importlib.import_module("planemo.galaxy.ephemeris_sleep"), "DEFAULT_SLEEP_WAIT", 0.05)

    started = time.monotonic()
    with pytest.raises(Exception, match="controlled startup failure") as exc_info:
        with serve_module.serve_daemon(
            config._ctx,
            daemon=True,
            port=port,
            skip_venv=True,
            galaxy_startup_timeout=200,
        ):
            pass

    assert "exited with code 7" in str(exc_info.value)
    assert time.monotonic() - started < 5
    assert config._daemon_process.returncode == 7
    assert not os.path.exists(config.pid_file)


def test_multiprocessing_daemon_survives_successful_planemo_exit(tmp_path, monkeypatch):
    galaxy_root = _make_galaxy_root(tmp_path, "25.0.1")
    _write_executable(galaxy_root / "run.sh", MODERN_RUN_SH)
    port = network_util.get_free_port()
    config = _make_local_config(galaxy_root, port)
    config.env["TEST_PORT"] = str(port)
    monkeypatch.setattr(serve_module, "galaxy_config", lambda *args, **kwds: _configured_galaxy(config))
    monkeypatch.setattr(importlib.import_module("planemo.galaxy.ephemeris_sleep"), "DEFAULT_SLEEP_WAIT", 0.05)

    try:
        with serve_module.serve(
            config._ctx,
            daemon=True,
            port=port,
            skip_venv=True,
            galaxy_startup_timeout=200,
        ):
            pass

        assert config._daemon_control_fd is None
        assert config._daemon_process.poll() is None
        assert network_util.wait_net_service("127.0.0.1", port, timeout=0.1)
    finally:
        config.kill()


PARENT_PROCESS = r"""
import json
import os
import shlex
import sys
import time
from types import SimpleNamespace

from planemo.galaxy.config import LocalGalaxyConfig

galaxy_root, config_directory, port, ready_path = sys.argv[1:]
env = {
    "GALAXY_LOG": os.path.join(galaxy_root, "main.log"),
    "GALAXY_PID": os.path.join(galaxy_root, "main.pid"),
    "GRAVITY_STATE_DIR": os.path.join(config_directory, "gravity"),
    "TEST_PORT": port,
    "TEST_RECORD": os.path.join(galaxy_root, "record.json"),
}
os.makedirs(env["GRAVITY_STATE_DIR"])
config = LocalGalaxyConfig(
    SimpleNamespace(verbose=False),
    config_directory,
    env,
    None,
    int(port),
    "main",
    "test-key",
    [],
    galaxy_root,
    {},
)
command = f"cd {shlex.quote(galaxy_root)} && ./run.sh --server-name main"
process = config.start_daemon(command)
with open(ready_path, "w") as ready_file:
    ready_file.write(str(process.pid))
while True:
    time.sleep(60)
"""


def test_sigkill_of_planemo_cleans_multiprocessing_process_group(tmp_path):
    galaxy_root = _make_galaxy_root(tmp_path, "25.0.1")
    _write_executable(galaxy_root / "run.sh", MODERN_RUN_SH)
    config_directory = galaxy_root / "planemo-config"
    config_directory.mkdir()
    port = network_util.get_free_port()
    ready_path = tmp_path / "parent-ready"
    record_path = galaxy_root / "record.json"
    repository_root = os.path.dirname(os.path.dirname(__file__))
    environ = os.environ.copy()
    environ["PYTHONPATH"] = repository_root
    parent = subprocess.Popen(
        [sys.executable, "-c", PARENT_PROCESS, str(galaxy_root), str(config_directory), str(port), str(ready_path)],
        env=environ,
    )
    monitor_pid = None
    try:
        _wait_for_path(ready_path)
        monitor_pid = int(ready_path.read_text())
        record = _wait_for_json(record_path)
        assert _pid_exists(monitor_pid)
        assert _pid_exists(record["pid"])
        assert _pid_exists(record["child_pid"])

        os.kill(parent.pid, signal.SIGKILL)
        parent.wait()

        _wait_for_pid_exit(monitor_pid)
        _wait_for_pid_exit(record["pid"])
        _wait_for_pid_exit(record["child_pid"])
    finally:
        if parent.poll() is None:
            parent.kill()
            parent.wait()
        if monitor_pid is not None and _pid_exists(monitor_pid):
            os.killpg(monitor_pid, signal.SIGKILL)


def test_monitor_escalates_for_sigterm_ignoring_process_group(tmp_path):
    galaxy_root = _make_galaxy_root(tmp_path, "25.0.1")
    _write_executable(galaxy_root / "run.sh", STUBBORN_RUN_SH)
    config_directory = galaxy_root / "planemo-config"
    config_directory.mkdir()
    port = network_util.get_free_port()
    ready_path = tmp_path / "parent-ready"
    record_path = galaxy_root / "record.json"
    repository_root = os.path.dirname(os.path.dirname(__file__))
    environ = os.environ.copy()
    environ["PYTHONPATH"] = repository_root
    parent = subprocess.Popen(
        [sys.executable, "-c", PARENT_PROCESS, str(galaxy_root), str(config_directory), str(port), str(ready_path)],
        env=environ,
    )
    monitor_pid = None
    try:
        _wait_for_path(ready_path)
        monitor_pid = int(ready_path.read_text())
        record = _wait_for_json(record_path)

        os.kill(parent.pid, signal.SIGKILL)
        parent.wait()

        _wait_for_pid_exit(monitor_pid)
        _wait_for_pid_exit(record["pid"])
        _wait_for_pid_exit(record["child_pid"])
    finally:
        if parent.poll() is None:
            parent.kill()
            parent.wait()
        if monitor_pid is not None and _pid_exists(monitor_pid):
            os.killpg(monitor_pid, signal.SIGKILL)
