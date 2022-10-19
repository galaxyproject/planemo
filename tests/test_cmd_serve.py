import os
import signal
import subprocess
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


class UsesServeCommand:
    def _run(self, serve_args=[], serve_cmd="serve"):
        serve_cmd = self._serve_command_list(serve_args, serve_cmd)
        if run_verbosely():
            print("Running command for test [%s]" % serve_cmd)
        if "--daemon" not in serve_cmd:
            self._run_subprocess(serve_cmd)
        else:
            self._check_exit_code(serve_cmd)

    def _run_subprocess(self, serve_cmd):
        serve_cmd.insert(0, "planemo")
        stdout_file = tempfile.NamedTemporaryFile(mode="wb+", suffix="_planemo_stdout")
        popen = subprocess.Popen(
            serve_cmd, env=os.environ.copy(), stdout=stdout_file, stderr=stdout_file, preexec_fn=os.setsid
        )
        if popen.poll() is not None:
            stdout_file.seek(0)
            raise Exception(
                f"planemo serve command failed. exit_code: {popen.returncode}, output: {stdout_file.read().decode('utf-8')}"
            )

        def cleanup():
            pgrp = os.getpgid(popen.pid)
            os.killpg(pgrp, signal.SIGINT)
            stdout_file.seek(0)
            print(stdout_file.read().decode("utf-8"))
            popen.terminate()
            popen.kill()

        self._cleanup_hooks.append(cleanup)

    def _serve_command_list(self, serve_args=[], serve_cmd="serve"):
        test_cmd = ["--verbose"] if run_verbosely() else []
        test_cmd.extend(
            [
                serve_cmd,
                "--galaxy_root",
                self.galaxy_root,
                "--galaxy_branch",
                target_galaxy_branch(),
                "--no_dependency_resolution",
                "--port",
                str(self._port),
                self._serve_artifact,
            ]
        )
        test_cmd.extend(serve_args)
        return test_cmd

    def _launch_thread_and_wait(self, func, args=None, run_as_subprocess=False, **kwd):
        args = args or []
        future = launch_and_wait_for_galaxy(self._port, func, [args], run_as_subprocess=run_as_subprocess, **kwd)
        if future:
            self._cleanup_hooks.append(lambda: future.cancel())

    @property
    def _user_gi(self):
        admin_gi = api.gi(self._port)
        user_api_key = api.user_api_key(admin_gi)
        user_gi = api.gi(self._port, key=user_api_key)
        return user_gi


class ServeTestCase(CliTestCase, UsesServeCommand):
    @classmethod
    def setUpClass(cls):
        cls.galaxy_root = tempfile.mkdtemp()

    @classmethod
    def tearDownClass(cls):
        safe_rmtree(cls.galaxy_root)

    def setUp(self):
        super().setUp()
        self._port = network_util.get_free_port()
        self._pid_file = os.path.join(self._home, "test.pid")
        self._serve_artifact = os.path.join(TEST_REPOS_DIR, "single_tool", "cat.xml")

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    @mark.tests_galaxy_branch
    def test_serve(self):
        extra_args = [
            "--skip_client_build",
        ]
        self._launch_thread_and_wait(self._run, extra_args, run_as_subprocess=True)

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    @skip_if_environ("PLANEMO_SKIP_GALAXY_CLIENT_TESTS")
    @skip_unless_executable("python3")
    def test_serve_client_python3(self):
        extra_args = ["--galaxy_python_version", "3"]
        # Given the client build - give this more time.
        timeout_multiplier = 3
        self._launch_thread_and_wait(
            self._run, extra_args, timeout_multiplier=timeout_multiplier, run_as_subprocess=True
        )
        # Check that the client was correctly built
        url = "http://localhost:%d/static/dist/analysis.bundled.js" % int(self._port)
        r = requests.get(url)
        assert r.status_code == 200

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    @mark.tests_galaxy_branch
    def test_serve_daemon(self):
        extra_args = ["--daemon", "--skip_client_build", "--pid_file", self._pid_file]
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
            "--pid_file",
            self._pid_file,
            "--extra_tools",
            random_lines,
            "--extra_tools",
            cat,
        ]
        self._launch_thread_and_wait(self._run, extra_args)
        time.sleep(30)
        user_gi = self._user_gi
        assert len(user_gi.histories.get_histories(name=TEST_HISTORY_NAME)) == 0
        user_gi.histories.create_history(TEST_HISTORY_NAME)
        assert user_gi.tools.show_tool("random_lines1")
        workflows = user_gi.workflows.get_workflows()
        assert len(workflows) == 1
        workflow = workflows[0]
        workflow_id = workflow["id"]
        assert workflow_id

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    @mark.tests_galaxy_branch
    def test_serve_multiple_tool_data_tables(self):
        tool_data_table_path_1 = self._tool_data_table("planemo1")
        tool_data_table_path_2 = self._tool_data_table("planemo2")
        extra_args = [
            "--daemon",
            "--skip_client_build",
            "--pid_file",
            self._pid_file,
            "--tool_data_table",
            tool_data_table_path_1,
            "--tool_data_table",
            tool_data_table_path_2,
        ]
        self._launch_thread_and_wait(self._run, extra_args)
        admin_gi = api.gi(self._port)
        table_contents = admin_gi.tool_data.show_data_table("__dbkeys__")
        assert any("planemo1" in field[0] for field in table_contents["fields"]), table_contents
        assert any("planemo2" in field[0] for field in table_contents["fields"]), table_contents

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    @skip_if_environ("PLANEMO_SKIP_SHED_TESTS")
    @mark.tests_galaxy_branch
    def test_shed_serve(self):
        extra_args = ["--daemon", "--skip_client_build", "--pid_file", self._pid_file, "--shed_target", "toolshed"]
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
            "--pid_file",
            self._pid_file,
            "--profile",
            new_profile,
        ]
        serve_cmd = self._serve_command_list(extra_args)
        with cli_daemon_galaxy(self._runner, self._pid_file, self._port, serve_cmd):
            user_gi = self._user_gi
            assert len(user_gi.histories.get_histories(name=TEST_HISTORY_NAME)) == 0
            user_gi.histories.create_history(TEST_HISTORY_NAME)

        # TODO: Pretty sure this is getting killed, but we should verify.
        with cli_daemon_galaxy(self._runner, self._pid_file, self._port, serve_cmd):
            assert len(user_gi.histories.get_histories(name=TEST_HISTORY_NAME)) == 1

    def _tool_data_table(self, dbkey):
        with tempfile.NamedTemporaryFile(
            "w", suffix=".xml.test", delete=False
        ) as tool_data_table, tempfile.NamedTemporaryFile("w", suffix="bla.loc", delete=False) as loc_file:
            tool_data_table.write(
                f"""<tables>
    <table name="__dbkeys__" comment_char="#">
        <columns>value, name, len_path</columns>
        <file path="{loc_file.name}" />
    </table>
</tables>"""
            )
            loc_file.write(f"{dbkey}\t{dbkey}\t{dbkey}.len")
            tool_data_table.flush()
            loc_file.flush()
        self._cleanup_hooks.extend([lambda: os.remove(loc_file.name), lambda: os.remove(tool_data_table.name)])
        return tool_data_table.name

    def _run_shed(self, serve_args=[]):
        return self._run(serve_args=serve_args, serve_cmd="shed_serve")
