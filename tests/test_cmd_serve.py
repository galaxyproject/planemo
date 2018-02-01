import functools
import os
import threading
import time
import uuid

from planemo import network_util
from planemo.galaxy import api
from planemo.io import kill_pid_file
from .test_utils import (
    CliTestCase,
    PROJECT_TEMPLATES_DIR,
    skip_if_environ,
    skip_unless_environ,
    TEST_DATA_DIR,
    TEST_REPOS_DIR,
)

TEST_HISTORY_NAME = "Cool History 42"


class ServeTestCase(CliTestCase):

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    def test_serve(self):
        self._launch_thread_and_wait(self._run)

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    def test_serve_daemon(self):
        extra_args = ["--daemon", "--pid_file", self._pid_file]
        self._launch_thread_and_wait(self._run, extra_args)
        user_gi = self._user_gi
        assert len(user_gi.histories.get_histories(name=TEST_HISTORY_NAME)) == 0
        user_gi.histories.create_history(TEST_HISTORY_NAME)

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    def test_serve_workflow(self):
        random_lines = os.path.join(PROJECT_TEMPLATES_DIR, "demo", "randomlines.xml")
        cat = os.path.join(PROJECT_TEMPLATES_DIR, "demo", "cat.xml")
        self._serve_artifact = os.path.join(TEST_DATA_DIR, "wf1.gxwf.yml")
        extra_args = [
            "--daemon",
            "--pid_file", self._pid_file,
            "--extra_tools", random_lines,
            "--extra_tools", cat,
        ]
        self._launch_thread_and_wait(self._run, extra_args)
        time.sleep(40)
        user_gi = self._user_gi
        assert len(user_gi.histories.get_histories(name=TEST_HISTORY_NAME)) == 0
        user_gi.histories.create_history(TEST_HISTORY_NAME)
        assert user_gi.tools.get_tools(tool_id="random_lines1")
        assert len(user_gi.workflows.get_workflows()) == 1

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    def test_serve_profile(self):
        self._test_serve_profile()

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    @skip_unless_environ("PLANEMO_ENABLE_POSTGRES_TESTS")
    def test_serve_postgres_profile(self):
        self._test_serve_profile("--database_type", "postgres")

    def _test_serve_profile(self, *db_options):
        new_profile = "planemo_test_profile_%s" % uuid.uuid4()
        extra_args = [
            "--daemon",
            "--pid_file", self._pid_file,
            "--profile", new_profile,
        ]
        extra_args.extend(db_options)
        self._launch_thread_and_wait(self._run, extra_args)
        user_gi = self._user_gi
        assert len(user_gi.histories.get_histories(name=TEST_HISTORY_NAME)) == 0
        user_gi.histories.create_history(TEST_HISTORY_NAME)
        kill_pid_file(self._pid_file)

        self._launch_thread_and_wait(self._run, extra_args)
        assert len(user_gi.histories.get_histories(name=TEST_HISTORY_NAME)) == 1

    def setUp(self):
        super(ServeTestCase, self).setUp()
        self._port = network_util.get_free_port()
        self._pid_file = os.path.join(self._home, "test.pid")
        self._serve_artifact = os.path.join(TEST_REPOS_DIR, "single_tool", "cat.xml")

    @property
    def _user_gi(self):
        admin_gi = api.gi(self._port)
        user_api_key = api.user_api_key(admin_gi)
        user_gi = api.gi(self._port, user_api_key)
        return user_gi

    def _launch_thread_and_wait(self, func, args=[]):
        target = functools.partial(func, args)
        port = self._port
        t = threading.Thread(target=target)
        t.daemon = True
        t.start()
        time.sleep(10)
        assert network_util.wait_net_service("127.0.0.1", port, timeout=600)
        time.sleep(1)
        assert network_util.wait_net_service("127.0.0.1", port, timeout=600)
        time.sleep(1)
        assert network_util.wait_net_service("127.0.0.1", port, timeout=600)

    def _run(self, serve_args=[]):
        test_cmd = [
            "serve",
            "--install_galaxy",
            "--port",
            str(self._port),
            self._serve_artifact,
        ]
        test_cmd.extend(serve_args)
        self._check_exit_code(test_cmd)
