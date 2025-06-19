"""The module contains a class to test the ``invocation_export`` command."""

import json
import os
import tempfile
import time
import zipfile

import pytest

from planemo import network_util
from .test_cmd_serve import UsesServeCommand
from .test_utils import (
    CliTestCase,
    mark,
    PROJECT_TEMPLATES_DIR,
    safe_rmtree,
    skip_if_environ,
    target_galaxy_branch,
    TEST_DATA_DIR,
)

TEST_HISTORY_NAME = "Cool History 42"
SERVE_TEST_VERBOSE = True


class GalaxyInvocationExportTestCase(CliTestCase, UsesServeCommand):
    @classmethod
    def setUpClass(cls):
        cls.galaxy_root = tempfile.mkdtemp()
        print(cls)

    @classmethod
    def tearDownClass(cls):
        safe_rmtree(cls.galaxy_root)

    def setUp(self):
        super().setUp()
        self._port = network_util.get_free_port()
        self._pid_file = os.path.join(self._home, "test.pid")

    @pytest.mark.skipif(
        target_galaxy_branch() == "release_22.05",
        reason="Skipping test on Galaxy 22.05, does not support rocrate.zip export.",
    )
    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    @mark.tests_galaxy_branch
    def test_export_invocation(self):
        with self._isolate() as f:
            # Create Galaxy instance
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
            workflows = user_gi.workflows.get_workflows()
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
                "--galaxy_user_key",
                user_gi.key,
                remote_uri,
            ]
            self._check_exit_code(test_command, exit_code=0)
            output_json_path = os.path.join(f, "tool_test_output.json")
            with open(output_json_path) as j_output:
                output = json.load(j_output)

            test_index = 1
            invocation_id = output["tests"][test_index]["data"]["invocation_details"]["details"]["invocation_id"]

            user_gi.invocations.wait_for_invocation(invocation_id)
            output_rocrate_path = os.path.join(f, "temp.rocrate.zip")
            export_invocation_cmd = [
                "invocation_export",
                "--galaxy_url",
                f"http://localhost:{self._port}",
                "--galaxy_user_key",
                user_gi.key,
                "--output",
                output_rocrate_path,
                invocation_id,
            ]
            self._check_exit_code(export_invocation_cmd, exit_code=0)

            assert zipfile.is_zipfile(output_rocrate_path)
