"""Keep a foreground Galaxy process tied to Planemo until detached."""

import os
import select
import signal
import subprocess
import sys
import time

DETACH = b"D"
POLL_INTERVAL = 0.1
TERMINATION_TIMEOUT = 0.5


def _process_group_exists(process_group):
    try:
        os.killpg(process_group, 0)
    except ProcessLookupError:
        return False
    return True


def _wait_for_process_group(process_group):
    deadline = time.monotonic() + TERMINATION_TIMEOUT
    while time.monotonic() < deadline:
        if not _process_group_exists(process_group):
            return True
        time.sleep(POLL_INTERVAL)
    return not _process_group_exists(process_group)


def _terminate_process_group(child):
    process_group = child.pid
    for process_signal in (signal.SIGTERM, signal.SIGKILL):
        try:
            os.killpg(process_group, process_signal)
        except ProcessLookupError:
            break
        if _wait_for_process_group(process_group):
            break
    child.wait()
    return child.returncode


def main():
    control_fd = int(sys.argv[1])
    command = sys.argv[2]
    child = subprocess.Popen(command, shell=True, start_new_session=True)
    termination_requested = False

    def request_termination(signum, frame):
        nonlocal termination_requested
        termination_requested = True

    signal.signal(signal.SIGTERM, request_termination)
    signal.signal(signal.SIGINT, request_termination)

    while child.poll() is None:
        if termination_requested:
            return _terminate_process_group(child)
        readable, _, _ = select.select([control_fd], [], [], POLL_INTERVAL)
        if not readable:
            continue
        message = os.read(control_fd, 1)
        if message == DETACH:
            os.close(control_fd)
            while child.poll() is None and not termination_requested:
                time.sleep(POLL_INTERVAL)
            if termination_requested:
                return _terminate_process_group(child)
            return child.returncode

        # EOF means Planemo exited without explicitly detaching the daemon.
        return _terminate_process_group(child)

    return child.returncode


if __name__ == "__main__":
    sys.exit(main())
