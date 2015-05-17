import functools
import os
import time
import threading

from .test_utils import (
    CliTestCase,
    skip_if_environ,
    TEST_REPOS_DIR,
)
from . import network_util


class ServeTestCase(CliTestCase):

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    def test_serve(self):
        port = network_util.get_free_port()
        serve = functools.partial(self._run, port)
        t = threading.Thread(target=serve)
        t.daemon = True
        t.start()
        time.sleep(15)
        assert network_util.wait_net_service("127.0.0.1", port)

    def _run(self, port):
        cat_path = os.path.join(TEST_REPOS_DIR, "single_tool", "cat.xml")
        test_cmd = [
            "serve",
            "--install_galaxy",
            "--port",
            str(port),
            cat_path,
        ]
        self._check_exit_code(test_cmd)
