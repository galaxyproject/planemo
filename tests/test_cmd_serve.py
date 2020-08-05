import os
import tempfile
import time
import uuid

import requests

from planemo import network_util
from planemo.galaxy import api
from planemo.io import kill_pid_file
from .test_utils import (
    cli_daemon_galaxy,
    CliTestCase,
    launch_and_wait_for_galaxy,
    mark,
    PROJECT_TEMPLATES_DIR,
    run_verbosely,
    safe_rmtree,
    skip_if_environ,
    skip_unless_environ,
    skip_unless_executable,
    target_galaxy_branch,
    TEST_DATA_DIR,
    TEST_REPOS_DIR,
)

TEST_HISTORY_NAME = "Cool History 42"
SERVE_TEST_VERBOSE = True


class ServeTestCase(CliTestCase):

    @classmethod
    def setUpClass(cls):
        cls.galaxy_root = tempfile.mkdtemp()

    @classmethod
    def tearDownClass(cls):
        safe_rmtree(cls.galaxy_root)

    def setUp(self):
        super(ServeTestCase, self).setUp()
        self._port = network_util.get_free_port()
        self._pid_file = os.path.join(self._home, "test.pid")
        self._serve_artifact = os.path.join(TEST_REPOS_DIR, "single_tool", "cat.xml")

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    @mark.tests_galaxy_branch
    def test_serve(self):
        extra_args = [
            "--skip_client_build",
        ]
        self._launch_thread_and_wait(self._run, extra_args)

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    @skip_if_environ("PLANEMO_SKIP_GALAXY_CLIENT_TESTS")
    @skip_unless_executable("python3")
    def test_serve_client_python3(self):
        extra_args = [
            "--galaxy_python_version", "3"
        ]
        # Given the client build - give this more time.
        timeout_multiplier = 3
        self._launch_thread_and_wait(self._run, extra_args, timeout_multiplier=timeout_multiplier)
        # Check that the client was correctly built
        url = "http://localhost:%d/static/dist/analysis.bundled.js" % int(self._port)
        r = requests.get(url)
        assert r.status_code == 200

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    @mark.tests_galaxy_branch
    def test_serve_daemon(self):
        extra_args = [
            "--daemon",
            "--skip_client_build",
            "--pid_file", self._pid_file]
        self._launch_thread_and_wait(self._run, extra_args)
        user_gi = self._user_gi
        assert len(user_gi.histories.get_histories(name=TEST_HISTORY_NAME)) == 0
        user_gi.histories.create_history(TEST_HISTORY_NAME)
        kill_pid_file(self._pid_file)

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    @mark.tests_galaxy_branch
    def test_serve_workflow(self):
        random_lines = os.path.join(PROJECT_TEMPLATES_DIR, "demo", "randomlines.xml")
        cat = os.path.join(PROJECT_TEMPLATES_DIR, "demo", "cat.xml")
        self._serve_artifact = os.path.join(TEST_DATA_DIR, "wf1.gxwf.yml")
        extra_args = [
            "--daemon",
            "--skip_client_build",
            "--pid_file", self._pid_file,
            "--extra_tools", random_lines,
            "--extra_tools", cat,
        ]
        self._launch_thread_and_wait(self._run, extra_args)
        time.sleep(30)
        user_gi = self._user_gi
        assert len(user_gi.histories.get_histories(name=TEST_HISTORY_NAME)) == 0
        user_gi.histories.create_history(TEST_HISTORY_NAME)
        assert user_gi.tools.get_tools(tool_id="random_lines1")
        assert len(user_gi.workflows.get_workflows()) == 1

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    @skip_if_environ("PLANEMO_SKIP_SHED_TESTS")
    @mark.tests_galaxy_branch
    def test_shed_serve(self):
        extra_args = [
            "--daemon",
            "--skip_client_build",
            "--pid_file", self._pid_file,
            "--shed_target", "toolshed"]
        fastqc_path = os.path.join(TEST_REPOS_DIR, "fastqc")
        self._serve_artifact = fastqc_path
        self._launch_thread_and_wait(self._run_shed, extra_args)
        user_gi = self._user_gi
        found = False
        tool_ids = None
        for i in range(30):
            tool_ids = [t["id"] for t in user_gi.tools.get_tools()]
            if any(_.startswith("toolshed.g2.bx.psu.edu/repos/devteam/fastqc/fastqc/") for _ in tool_ids):
                found = True
                break
            time.sleep(5)

        assert found, "Failed to find fastqc id in %s" % tool_ids
        kill_pid_file(self._pid_file)

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    def test_serve_profile(self):
        self._test_serve_profile()

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    @skip_unless_environ("PLANEMO_ENABLE_POSTGRES_TESTS")
    @skip_unless_executable("psql")
    def test_serve_postgres_profile(self):
        self._test_serve_profile("--database_type", "postgres")

    def _test_serve_profile(self, *db_options):
        new_profile = "planemo_test_profile_%s" % uuid.uuid4()
        extra_args = [
            "--daemon",
            "--skip_client_build",
            "--pid_file", self._pid_file,
            "--profile", new_profile,
        ]
        serve_cmd = self._serve_command_list(extra_args)
        with cli_daemon_galaxy(self._runner, self._pid_file, self._port, serve_cmd):
            user_gi = self._user_gi
            assert len(user_gi.histories.get_histories(name=TEST_HISTORY_NAME)) == 0
            user_gi.histories.create_history(TEST_HISTORY_NAME)

        # TODO: Pretty sure this is getting killed, but we should verify.
        with cli_daemon_galaxy(self._runner, self._pid_file, self._port, serve_cmd):
            assert len(user_gi.histories.get_histories(name=TEST_HISTORY_NAME)) == 1

    @property
    def _user_gi(self):
        admin_gi = api.gi(self._port)
        user_api_key = api.user_api_key(admin_gi)
        user_gi = api.gi(self._port, key=user_api_key)
        return user_gi

    def _launch_thread_and_wait(self, func, args=[], **kwd):
        future = launch_and_wait_for_galaxy(self._port, func, [args], **kwd)
        if future is not None:
            self._futures.append(future)

    def _run_shed(self, serve_args=[]):
        return self._run(serve_args=serve_args, serve_cmd="shed_serve")

    def _run(self, serve_args=[], serve_cmd="serve"):
        serve_cmd = self._serve_command_list(serve_args, serve_cmd)
        if run_verbosely():
            print("Running command for test [%s]" % serve_cmd)
        self._check_exit_code(serve_cmd)

    def _serve_command_list(self, serve_args=[], serve_cmd="serve"):
        test_cmd = ["--verbose"] if run_verbosely() else []
        test_cmd.extend([
            serve_cmd,
            "--galaxy_root", self.galaxy_root,
            "--galaxy_branch", target_galaxy_branch(),
            "--no_dependency_resolution",
            "--port", str(self._port),
            self._serve_artifact,
        ])
        test_cmd.extend(serve_args)
        return test_cmd
