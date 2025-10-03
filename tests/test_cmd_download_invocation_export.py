"""Module contains :class:`CmdTestTestCase` - integration tests for the ``test`` command."""

import os
import tempfile
import time

from planemo import network_util
from .test_cmd_serve import UsesServeCommand
from .test_utils import (
    CliTestCase,
    PROJECT_TEMPLATES_DIR,
    safe_rmtree,
    skip_if_environ,
    TEST_DATA_DIR,
)


class CmdTestTestCase(CliTestCase, UsesServeCommand):
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
    def test_download_run_output(self):
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

            # Run workflow and download the output
            wf = os.path.join(TEST_DATA_DIR, "wf2.ga")
            job_path = os.path.join(TEST_DATA_DIR, "wf2-job.yml")
            output_path_org = os.path.join(f, "output_org")
            if not os.path.exists(output_path_org):
                os.makedirs(output_path_org)
            output_path_download = os.path.join(f, "output_download")
            if not os.path.exists(output_path_download):
                os.makedirs(output_path_download)
            run_workflow_cmd = [
                "--verbose",
                "run",
                "--engine",
                "external_galaxy",
                "--galaxy_url",
                f"http://localhost:{self._port}",
                "--galaxy_user_key",
                user_gi.key,
                "--download_outputs",
                "--output_directory",
                output_path_org,
                wf,
                job_path,
            ]

            self._check_exit_code(run_workflow_cmd)

            # Get invocation and download generated output
            invocations = user_gi.invocations.get_invocations()
            download_outputs_cmd = [
                "--verbose",
                "invocation_download",
                "--galaxy_url",
                f"http://localhost:{self._port}",
                "--galaxy_user_key",
                user_gi.key,
                "--output_directory",
                output_path_download,
                invocations[0]["id"],
            ]
            self._check_exit_code(download_outputs_cmd)

            output_download_run_files = os.listdir(output_path_download)
            output_run_files = os.listdir(output_path_org)

            # Compare downloading result at runtime with download the output
            # after the run has completed using a separate command.
            assert len(output_download_run_files) == len(output_run_files)
            for i in range(len(output_download_run_files)):
                assert output_download_run_files[i] == output_run_files[i]
