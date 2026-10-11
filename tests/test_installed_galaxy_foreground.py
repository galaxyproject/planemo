"""Exercise foreground ownership with real processes without installing Galaxy."""

import json
import os
import signal
import subprocess
import sys
import time

import pytest
from galaxy.util.commands import argv_to_str

from planemo.galaxy.config import installed_galaxy_config
from planemo.io import (
    process_group_exists,
    terminate_process_group,
    TERMINATION_TIMEOUT_ENVIRON_KEY,
)
from .test_utils import create_test_context

SERVICE_SCRIPT = """
import json
import os
import signal
import subprocess
import sys
import time

child = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(300)"])
def stop(signum, frame):
    child.wait(timeout=5)
    raise SystemExit(0)
signal.signal(signal.SIGTERM, stop)
with open(sys.argv[1] + ".tmp", "w") as record:
    json.dump({"group": os.getpgrp(), "child": child.pid}, record)
os.replace(sys.argv[1] + ".tmp", sys.argv[1])
while True:
    time.sleep(1)
"""


@pytest.mark.parametrize("termination_signal", [signal.SIGTERM, signal.SIGINT, signal.SIGKILL])
def test_foreground_parent_exit_stops_service_group(tmp_path, termination_signal):
    service = tmp_path / "service.py"
    service.write_text(SERVICE_SCRIPT)
    record = tmp_path / "service.json"
    command = argv_to_str([sys.executable, str(service), str(record)])
    driver = (
        "from planemo.galaxy.config import installed_galaxy_config\n"
        "from tests.test_utils import create_test_context\n"
        "with installed_galaxy_config(create_test_context(), [], port=8765) as config:\n"
        f"    config.run_foreground({command!r})\n"
    )
    environment = os.environ.copy()
    environment[TERMINATION_TIMEOUT_ENVIRON_KEY] = "2"
    group = None
    with open(tmp_path / "foreground.log", "w+") as log:
        process = subprocess.Popen(
            [sys.executable, "-c", driver], env=environment, start_new_session=True, stdout=log, stderr=log
        )
        try:
            deadline = time.monotonic() + 30
            while not record.exists() and process.poll() is None and time.monotonic() < deadline:
                time.sleep(0.05)
            log.seek(0)
            assert record.exists(), log.read()
            group = json.loads(record.read_text())["group"]
            assert group != process.pid
            os.killpg(process.pid, termination_signal)
            process.wait(timeout=15)
            assert process.returncode != 0
            deadline = time.monotonic() + 15
            while process_group_exists(group) and time.monotonic() < deadline:
                time.sleep(0.05)
            assert not process_group_exists(group)
        finally:
            if group is not None:
                terminate_process_group(group, timeout=2)
            terminate_process_group(process.pid, timeout=2, reap=process.poll)


@pytest.mark.parametrize("exit_code", [0, 7])
def test_foreground_preserves_exit_code_and_terminal_output(capfd, exit_code):
    command = argv_to_str(
        [
            sys.executable,
            "-c",
            f"import sys; print('service stdout'); print('service stderr', file=sys.stderr); sys.exit({exit_code})",
        ]
    )
    with installed_galaxy_config(create_test_context(), [], port=8765) as config:
        assert config.run_foreground(command) == exit_code
        assert config._daemon_control_fd is None
        assert config._daemon_process.poll() == exit_code
    stdout, stderr = capfd.readouterr()
    assert "service stdout" in stdout
    assert "service stderr" in stderr
