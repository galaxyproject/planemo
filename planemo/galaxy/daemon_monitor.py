"""Keep a foreground Galaxy process tied to Planemo until detached.

Planemo runs this as ``python -m planemo.galaxy.daemon_monitor``, which imports
the ``planemo.galaxy`` package first - so all of Planemo is already loaded here
and sharing its process-group helpers costs nothing.
"""

import os
import select
import signal
import subprocess
import sys
import time

from planemo.io import terminate_process_group

DETACH = b"D"
POLL_INTERVAL = 0.1


def _terminate_process_group(child):
    """Stop Galaxy's process group and reap its leader."""
    # start_new_session=True made the child a process group leader, so its PID
    # doubles as a group ID that stays valid even after the leader is reaped.
    if terminate_process_group(child.pid, reap=child.poll):
        return child.wait()
    return 1


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
