from collections import namedtuple
import contextlib
import shutil
import socket
from tempfile import mkdtemp
import threading
from requests import post
from time import time as now

from werkzeug.serving import run_simple

from .shed_app import (
    app,
    InMemoryShedDataModel,
)

DEFAULT_OP_TIMEOUT = 2


def mock_model(directory):
    return InMemoryShedDataModel(
        directory
    ).add_category(
        "c1", "Text Manipulation"
    ).add_category(
        "c2", "Sequence Analysis"
    ).add_repository(
        "r1",
        name="test_repo_1",
        owner="iuc",
    )


def setup_mock_shed():
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    sock.bind(('localhost', 0))
    port = sock.getsockname()[1]
    sock.close()
    directory = mkdtemp()

    def run():
        app.debug = True
        app.config["model"] = mock_model(directory)
        run_simple(
            'localhost',
            port,
            app,
            use_reloader=False,
            use_debugger=True
        )

    t = threading.Thread(target=run)
    t.start()
    _wait_net_service("localhost", port, timeout=DEFAULT_OP_TIMEOUT)
    return MockShed("http://localhost:%d" % port, directory, t)


@contextlib.contextmanager
def mock_shed():
    mock_shed_obj = None
    try:
        mock_shed_obj = setup_mock_shed()
        yield mock_shed_obj
    finally:
        if mock_shed_obj is not None:
            mock_shed_obj.shutdown()


# code.activestate.com/recipes/576655-wait-for-network-service-to-appear
def _wait_net_service(server, port, timeout=None):
    """ Wait for network service to appear
        @param timeout: in seconds, if None or 0 wait forever
        @return: True of False, if timeout is None may return only True or
                 throw unhandled network exception
    """
    s = socket.socket()
    if timeout:
        end = now() + timeout

    while True:
        try:
            if timeout:
                next_timeout = end - now()
                if next_timeout < 0:
                    return False
                else:
                    s.settimeout(next_timeout)

            s.connect((server, port))

        except socket.timeout:
            # this exception occurs only if timeout is set
            if timeout:
                return False

        except socket.error:
            pass
        else:
            s.close()
            return True


def _shutdown(self):
    post("%s/shutdown" % self.url)
    self.thread.join(DEFAULT_OP_TIMEOUT)
    shutil.rmtree(self.directory)

MockShed = namedtuple("MockShed", ["url", "directory", "thread"])
MockShed.shutdown = _shutdown

__all__ = ["setup_mock_shed", "mock_shed"]
