import functools
import os
import time
import threading

from planemo.galaxy import api

from .test_utils import (
    CliTestCase,
    skip_if_environ,
    TEST_REPOS_DIR,
)
from . import network_util

TEST_HISTORY_NAME = "Cool History 42"


class ServeTestCase(CliTestCase):

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    def test_serve(self):
        port = network_util.get_free_port()
        serve = functools.partial(self._run, port)
        self._launch_thread_and_wait(serve, port)

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    def test_serve_daemon(self):
        port = network_util.get_free_port()
        pid_file = os.path.join(self._home, "test.pid")
        extra_args = ["--daemon", "--pid_file", pid_file]
        serve = functools.partial(self._run, port, extra_args)
        self._launch_thread_and_wait(serve, port)
        time.sleep(.1)
        assert network_util.wait_net_service("127.0.0.1", port)
        admin_gi = api.gi(port)
        user_api_key = api.user_api_key(admin_gi)
        user_gi = api.gi(port, user_api_key)
        assert len(user_gi.histories.get_histories(name=TEST_HISTORY_NAME)) == 0
        user_gi.histories.create_history(TEST_HISTORY_NAME)

    def _launch_thread_and_wait(self, func, port):
        t = threading.Thread(target=func)
        t.daemon = True
        t.start()
        time.sleep(15)
        assert network_util.wait_net_service("127.0.0.1", port)

    def _run(self, port, serve_args=[]):
        cat_path = os.path.join(TEST_REPOS_DIR, "single_tool", "cat.xml")
        test_cmd = [
            "serve",
            "--install_galaxy",
            "--port",
            str(port),
            cat_path,
        ]
        test_cmd.extend(serve_args)
        self._check_exit_code(test_cmd)
