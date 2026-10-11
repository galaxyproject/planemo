"""Opt-in integration coverage for package-installed Galaxy under Gravity.

Run with ``PLANEMO_TEST_INSTALLED_GALAXY=1`` in an environment containing a
supported Galaxy package set.
"""

import json
import os
import shutil
import signal
import socket

import pytest

from planemo import network_util
from planemo.galaxy.ephemeris_sleep import sleep
from planemo.io import process_group_exists
from .test_utils import (
    CliTestCase,
    planemo_subprocess,
    PROJECT_TEMPLATES_DIR,
    skip_unless_environ,
    TEST_DATA_DIR,
    TEST_TOOLS_DIR,
)

SUBPROCESS_STARTUP_TIMEOUT = 120
SUBPROCESS_TEST_TIMEOUT = 300
SUBPROCESS_TOOL_SHED_TIMEOUT = 600
SUBPROCESS_EXIT_TIMEOUT = 45


@pytest.mark.installed_galaxy
class InstalledGalaxyIntegrationTestCase(CliTestCase):
    @skip_unless_environ("PLANEMO_TEST_INSTALLED_GALAXY")
    def test_cli_test_covers_xml_yaml_and_local_workflow(self):
        xml_tool_path = os.path.join(PROJECT_TEMPLATES_DIR, "demo", "cat.xml")
        yaml_tool_path = os.path.join(TEST_TOOLS_DIR, "installed_echo.yml")
        workflow_path = os.path.join(TEST_DATA_DIR, "wf2.ga")

        with self._isolate() as test_directory:
            report_path = os.path.join(test_directory, "installed-report.json")
            process_log = os.path.join(test_directory, "planemo-test.log")
            command = [
                "test",
                "--engine",
                "installed_galaxy",
                "--no_dependency_resolution",
                "--extra_tools",
                xml_tool_path,
                "--test_output_json",
                report_path,
                xml_tool_path,
                yaml_tool_path,
                workflow_path,
            ]
            with planemo_subprocess(command, test_directory, process_log) as managed_process:
                process = managed_process.process
                return_code = process.wait(timeout=SUBPROCESS_TEST_TIMEOUT)
                process_output = managed_process.read_output()
                assert return_code == 0, process_output
                assert not process_group_exists(process.pid), process_output

            with open(report_path) as report_fh:
                report = json.load(report_fh)

        assert report["summary"] == {
            "num_errors": 0,
            "num_failures": 0,
            "num_skips": 0,
            "num_tests": 3,
        }
        assert all(test["data"]["status"] == "success" for test in report["tests"])
        tool_jobs = [test["data"]["job"] for test in report["tests"] if test["data"].get("job")]
        assert {job["tool_id"] for job in tool_jobs} == {"cat", "installed_echo"}
        assert all(job["state"] == "ok" for job in tool_jobs)
        # Pydantic's report model serializes optional fields as null, so select
        # by the authoritative runnable type rather than key presence.
        workflow_results = [test for test in report["tests"] if test["test_type"] == "galaxy_workflow"]
        assert len(workflow_results) == 1
        assert all(result["data"]["invocation_details"] for result in workflow_results)
        assert {result["data"]["invocation_details"]["details"]["invocation_state"] for result in workflow_results} <= {
            "scheduled",
            "completed",
        }

    @skip_unless_environ("PLANEMO_TEST_INSTALLED_GALAXY")
    @pytest.mark.installed_galaxy_toolshed
    def test_cli_test_installs_tool_shed_workflow_dependencies(self):
        shed_workflow_source = os.path.join(
            TEST_DATA_DIR,
            "wf_repos",
            "basic_wf_iwc_invalid_version",
            "Super-simple-workflow.ga",
        )

        with self._isolate() as test_directory:
            # The source is a workflow-lint fixture with an intentionally
            # invalid tool version and a content-specific test. Keep this
            # acceptance test focused on repairing/installing the Tool Shed
            # dependency and executing it, without changing that fixture's
            # separate contract.
            shed_workflow_path = os.path.join(test_directory, "installed-shed-workflow.ga")
            shutil.copyfile(shed_workflow_source, shed_workflow_path)
            shed_tests_path = shed_workflow_path.replace(".ga", "-tests.yml")
            with open(shed_tests_path, "w") as shed_tests:
                shed_tests.write("""- doc: Exercise an installed Tool Shed workflow
  job:
    n_rows: 5
  outputs:
    outfile:
      asserts:
        has_n_lines:
          n: 5
""")
            report_path = os.path.join(test_directory, "installed-shed-report.json")
            process_log = os.path.join(test_directory, "planemo-shed-test.log")
            command = [
                "test",
                "--engine",
                "installed_galaxy",
                "--no_dependency_resolution",
                "--test_output_json",
                report_path,
                shed_workflow_path,
            ]
            with planemo_subprocess(command, test_directory, process_log) as managed_process:
                process = managed_process.process
                return_code = process.wait(timeout=SUBPROCESS_TOOL_SHED_TIMEOUT)
                process_output = managed_process.read_output()
                assert return_code == 0, process_output
                assert not process_group_exists(process.pid), process_output

            with open(report_path) as report_fh:
                report = json.load(report_fh)

        assert report["summary"] == {
            "num_errors": 0,
            "num_failures": 0,
            "num_skips": 0,
            "num_tests": 1,
        }
        assert report["tests"][0]["test_type"] == "galaxy_workflow"
        assert report["tests"][0]["data"]["status"] == "success"

    @skip_unless_environ("PLANEMO_TEST_INSTALLED_GALAXY")
    def test_repeated_runs_use_fresh_gravity_process_groups(self):
        yaml_tool_path = os.path.join(TEST_TOOLS_DIR, "installed_echo.yml")

        with self._isolate() as test_directory:
            for invocation in range(2):
                expected_message = f"hello from planemo run {invocation}"
                job_path = os.path.join(test_directory, f"installed-job-{invocation}.json")
                output_directory = os.path.join(test_directory, f"outputs-{invocation}")
                output_json = os.path.join(test_directory, f"run-outputs-{invocation}.json")
                process_log = os.path.join(test_directory, f"planemo-run-{invocation}.log")
                with open(job_path, "w") as job_fh:
                    json.dump({"message": expected_message}, job_fh)

                command = [
                    "run",
                    "--engine",
                    "installed_galaxy",
                    "--no_dependency_resolution",
                    "--download_outputs",
                    "--output_directory",
                    output_directory,
                    "--output_json",
                    output_json,
                    yaml_tool_path,
                    job_path,
                ]
                with planemo_subprocess(command, test_directory, process_log) as managed_process:
                    process = managed_process.process
                    return_code = process.wait(timeout=SUBPROCESS_STARTUP_TIMEOUT)
                    process_output = managed_process.read_output()
                    assert return_code == 0, process_output
                    assert not process_group_exists(process.pid), process_output

                with open(output_json) as output_fh:
                    run_outputs = json.load(output_fh)
                assert "output" in run_outputs
                output_files = [path for path in os.scandir(output_directory) if path.is_file()]
                assert len(output_files) == 1
                with open(output_files[0].path) as output_fh:
                    assert output_fh.read().strip() == expected_message

    @skip_unless_environ("PLANEMO_TEST_INSTALLED_GALAXY")
    def test_foreground_serve_exits_cleanly_on_sigint(self):
        yaml_tool_path = os.path.join(TEST_TOOLS_DIR, "installed_echo.yml")

        with self._isolate() as test_directory:
            port = network_util.get_free_port()
            process_log = os.path.join(test_directory, "planemo-serve.log")
            command = [
                "serve",
                "--engine",
                "installed_galaxy",
                "--no_dependency_resolution",
                "--port",
                str(port),
                yaml_tool_path,
            ]
            with planemo_subprocess(command, test_directory, process_log) as managed_process:
                process = managed_process.process
                galaxy_url = f"http://127.0.0.1:{port}"
                if not sleep(galaxy_url, timeout=SUBPROCESS_STARTUP_TIMEOUT):
                    raise AssertionError(f"Installed Galaxy did not become ready.\n{managed_process.read_output()}")

                # Signal Planemo's foreground group as a terminal does. Gravity
                # runs in a child group, so Planemo must propagate teardown and
                # wait for every service before returning to the shell.
                os.killpg(process.pid, signal.SIGINT)
                return_code = process.wait(timeout=SUBPROCESS_EXIT_TIMEOUT)
                process_output = managed_process.read_output()
                assert return_code == 1, process_output
                assert "Aborted!" in process_output
                assert not process_group_exists(process.pid), process_output
                with pytest.raises(OSError):
                    socket.create_connection(("127.0.0.1", port), timeout=1)
