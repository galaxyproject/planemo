import json
import os
import tempfile
import time

from planemo import network_util
from planemo.galaxy import api
from .test_cmd_serve import UsesServeCommand
from .test_utils import (
    CliTestCase,
    mark,
    PROJECT_TEMPLATES_DIR,
    safe_rmtree,
    skip_if_environ,
    TEST_DATA_DIR,
)

TEST_HISTORY_NAME = "Cool History 42"
SERVE_TEST_VERBOSE = True


class CmdsWithWorkflowIdTestCase(CliTestCase, UsesServeCommand):
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

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    @mark.tests_galaxy_branch
    def test_serve_workflow(self):
        with self._isolate() as f:
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

            pseudo_path = os.path.join(TEST_DATA_DIR, "wf11-remote.gxwf.yml")
            remote_uri = f"gxid://workflows/{workflow_id}?runnable_path={pseudo_path}"
            test_command = [
                "test",
                "--engine",
                "external_galaxy",
                "--galaxy_url",
                f"http://localhost:{self._port}",
                "--galaxy_admin_key",
                api.DEFAULT_ADMIN_API_KEY,
                remote_uri,
            ]
            self._check_exit_code(test_command, exit_code=0)
            output_json_path = os.path.join(f, "tool_test_output.json")
            with open(output_json_path) as f:
                output = json.load(f)
            assert "tests" in output
