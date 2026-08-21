"""Test utilities from :module:`planemo.io`."""

import contextlib
import signal
import subprocess
import sys
import tempfile
import time

from planemo import io
from .test_utils import assert_equal


def test_io_capture():
    """Test :func:`planemo.io.conditionally_captured_io`."""
    with io.conditionally_captured_io(True, tee=False) as capture:
        io.warn("Problem...")
    assert_equal(capture[0]["data"], "Problem...")

    with io.conditionally_captured_io(True, tee=False) as capture:
        io.shell("echo 'Problem...'")
    assert_equal(capture[0]["data"], "echo 'Problem...'")
    assert_equal(capture[1]["data"], "Problem...")

    with io.conditionally_captured_io(True, tee=False) as capture:
        io.communicate("echo 'Problem...'")
    assert_equal(capture[0]["data"], "echo 'Problem...'")
    assert_equal(capture[1]["data"], "Problem...")

    with io.conditionally_captured_io(False, tee=False) as capture:
        io.communicate("echo 'Test...'")

    assert capture is None


def test_filter_paths():
    """Test :func:`planemo.io.filter_paths`."""
    test_cwd = "/a/b"

    def assert_filtered_is(paths, expected, **kwds):
        result = io.filter_paths(paths, cwd=test_cwd, **kwds)
        assert result == expected, f"paths [{result}] arent't expected [{expected}]"

    assert_filtered_is([], [], exclude=["/a"])
    assert_filtered_is(["/a/c"], [], exclude=["/a"])
    assert_filtered_is(["/b"], ["/b"], exclude=["/a"])
    assert_filtered_is(["/a/b/c"], [], exclude=["c"])
    with tempfile.NamedTemporaryFile(mode="w+") as tmp:
        tmp.write("#exclude c\n\nc\n")
        tmp.flush()
        assert_filtered_is(["/a/b/c", "/a/b/d"], ["/a/b/d"], exclude_from=[tmp.name])


SIGTERM_IGNORING_PROCESS = """
import signal
import sys
import time

signal.signal(signal.SIGTERM, signal.SIG_IGN)
with open(sys.argv[1], "w") as ready_file:
    ready_file.write("ready")
time.sleep(300)
"""


def _spawn_sigterm_ignoring_group(ready_path):
    """Start a leader that ignores SIGTERM, in a process group of its own.

    Waits for the handler to actually be installed - signalling before that
    point would kill the process outright and prove nothing.
    """
    process = subprocess.Popen(
        [sys.executable, "-c", SIGTERM_IGNORING_PROCESS, str(ready_path)],
        start_new_session=True,
    )
    for _ in range(200):
        if ready_path.exists():
            return process
        time.sleep(0.05)
    process.kill()
    process.wait()
    raise AssertionError("Process never began ignoring SIGTERM")


@contextlib.contextmanager
def _sigterm_ignoring_group(ready_path):
    process = _spawn_sigterm_ignoring_group(ready_path)
    try:
        yield process
    finally:
        if process.poll() is None:
            process.kill()
            process.wait()


def test_terminate_process_group_escalates_to_sigkill(tmp_path):
    """A group that ignores SIGTERM is still killed."""
    with _sigterm_ignoring_group(tmp_path / "ready") as process:
        assert io.terminate_process_group(process.pid, timeout=0.2, reap=process.poll)
        assert not io.process_group_exists(process.pid)
        assert process.returncode == -signal.SIGKILL


def test_terminate_process_group_does_not_escalate_for_a_willing_process(tmp_path):
    """A group that honours SIGTERM is never SIGKILLed."""
    process = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(300)"], start_new_session=True)
    try:
        assert io.terminate_process_group(process.pid, timeout=5, reap=process.poll)
        assert process.returncode == -signal.SIGTERM
    finally:
        if process.poll() is None:
            process.kill()
            process.wait()


def test_kill_posix_escalates_to_sigkill(tmp_path):
    """:func:`planemo.io.kill_posix` escalates rather than giving up on SIGTERM."""
    with _sigterm_ignoring_group(tmp_path / "ready") as process:
        io.kill_posix(process.pid)
        assert process.wait(timeout=5) == -signal.SIGKILL


def test_kill_posix_tolerates_missing_process():
    """Killing an already-reaped pid is a no-op rather than an error."""
    process = subprocess.Popen([sys.executable, "-c", ""], start_new_session=True)
    process.wait()
    io.kill_posix(process.pid)
