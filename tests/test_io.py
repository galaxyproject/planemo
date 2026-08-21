"""Test utilities from :module:`planemo.io`."""

import signal
import subprocess
import sys
import tempfile

from planemo import io
from .test_utils import (
    assert_equal,
    sigterm_ignoring_group,
)


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


def test_terminate_process_group_escalates_to_sigkill(tmp_path):
    """A group that ignores SIGTERM is still killed."""
    with sigterm_ignoring_group(tmp_path / "ready") as process:
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


def test_kill_posix_escalates_to_sigkill(tmp_path, monkeypatch):
    """:func:`planemo.io.kill_posix` escalates rather than giving up on SIGTERM."""
    monkeypatch.setenv(io.TERMINATION_TIMEOUT_ENVIRON_KEY, "0.5")
    with sigterm_ignoring_group(tmp_path / "ready") as process:
        io.kill_posix(process.pid)
        assert process.wait(timeout=5) == -signal.SIGKILL


def test_kill_posix_tolerates_missing_process():
    """Killing an already-reaped pid is a no-op rather than an error."""
    process = subprocess.Popen([sys.executable, "-c", ""], start_new_session=True)
    process.wait()
    io.kill_posix(process.pid)


def test_termination_timeout_honours_the_environment(monkeypatch):
    """The grace period is tunable, and a bad value falls back rather than raising."""
    monkeypatch.delenv(io.TERMINATION_TIMEOUT_ENVIRON_KEY, raising=False)
    assert io.termination_timeout() == io.DEFAULT_TERMINATION_TIMEOUT
    monkeypatch.setenv(io.TERMINATION_TIMEOUT_ENVIRON_KEY, "0.25")
    assert io.termination_timeout() == 0.25
    monkeypatch.setenv(io.TERMINATION_TIMEOUT_ENVIRON_KEY, "soon")
    assert io.termination_timeout() == io.DEFAULT_TERMINATION_TIMEOUT
