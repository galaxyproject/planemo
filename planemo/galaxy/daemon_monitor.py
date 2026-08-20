"""Keep a foreground Galaxy process tied to Planemo until detached."""

import os
import select
import signal
import subprocess
import sys

DETACH = b"D"
POLL_INTERVAL = 0.1


def main():
    control_fd = int(sys.argv[1])
    command = sys.argv[2]
    child = subprocess.Popen(command, shell=True)

    while child.poll() is None:
        readable, _, _ = select.select([control_fd], [], [], POLL_INTERVAL)
        if not readable:
            continue
        message = os.read(control_fd, 1)
        if message == DETACH:
            os.close(control_fd)
            return child.wait()

        # EOF means Planemo exited without explicitly detaching the daemon.
        # The monitor and Galaxy share a dedicated process group, so signaling
        # the group also cleans up every service launched by Gravity.
        os.killpg(os.getpgrp(), signal.SIGTERM)

    return child.returncode


if __name__ == "__main__":
    sys.exit(main())
