"""Module contains :class:`CmdTestTestCase` - integration tests for the ``test`` command."""
import json
import os

from .test_utils import (
    assert_exists,
    CliTestCase,
    PROJECT_TEMPLATES_DIR,
    skip_if_environ,
    skip_unless_module,
    TEST_DATA_DIR,
)


class CmdTestTestCase(CliTestCase):
    """Integration tests for the ``test`` command."""

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    def test_workflow_test_simple_yaml(self):
        """Test testing a simple YAML workflow with Galaxy."""
        with self._isolate():
            random_lines = os.path.join(PROJECT_TEMPLATES_DIR, "demo", "randomlines.xml")
            cat = os.path.join(PROJECT_TEMPLATES_DIR, "demo", "cat.xml")
            test_artifact = os.path.join(TEST_DATA_DIR, "wf1.gxwf.yml")
            test_command = [
                "--verbose",
                "test",
            ]
            test_command = self.append_profile_argument_if_needed(test_command)
            test_command += [
                "--no_dependency_resolution",
                "--galaxy_branch", "release_18.01",  # Much better workflow output detection than master for now (pre-release of 18.01)
                "--extra_tools", random_lines,
                "--extra_tools", cat,
                test_artifact,
            ]
            self._check_exit_code(test_command, exit_code=0)

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    def test_workflow_test_simple_ga(self):
        """Test testing a simple GA workflow with Galaxy."""
        with self._isolate():
            cat = os.path.join(PROJECT_TEMPLATES_DIR, "demo", "cat.xml")
            test_artifact = os.path.join(TEST_DATA_DIR, "wf2.ga")
            test_command = [
                "--verbose",
                "test"
            ]
            test_command = self.append_profile_argument_if_needed(test_command)
            test_command += [
                "--no_dependency_resolution",
                "--galaxy_branch", "release_18.01",  # Much better workflow output detection than master for now (pre-release of 18.01)
                "--extra_tools", cat,
                test_artifact,
            ]
            # try:
            self._check_exit_code(test_command, exit_code=0)
            # except Exception:
            #    with open(os.path.join(f, "tool_test_output.json"), "r") as o:
            #        print(o.read())
            #    raise

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    def test_workflow_test_distro_tool(self):
        """Test testing a simple GA workflow with Galaxy."""
        with self._isolate():
            test_artifact = os.path.join(TEST_DATA_DIR, "wf4-distro-tools.gxwf.yml")
            test_command = [
                "--verbose",
                "test"
            ]
            test_command = self.append_profile_argument_if_needed(test_command)
            test_command += [
                "--no_dependency_resolution",
                "--galaxy_branch", "release_18.01",  # Much better workflow output detection than master for now (pre-release of 18.01)
                test_artifact,
            ]
            # try:
            self._check_exit_code(test_command, exit_code=0)
            # except Exception:
            #    with open(os.path.join(f, "tool_test_output.json"), "r") as o:
            #        print(o.read())
            #    raise

    @skip_if_environ("PLANEMO_SKIP_CWLTOOL_TESTS")
    def test_cwltool_tool_test(self):
        """Test testing a CWL tool with cwltool."""
        with self._isolate() as f:
            test_artifact = os.path.join(TEST_DATA_DIR, "int_tool.cwl")
            test_command = [
                "test",
                test_artifact,
            ]
            self._check_exit_code(test_command, exit_code=0)
            assert_exists(os.path.join(f, "tool_test_output.json"))

    @skip_if_environ("PLANEMO_SKIP_CWLTOOL_TESTS")
    def test_cwltool_tool_url_inputs_test(self):
        """Test testing a CWL tool with cwltool."""
        with self._isolate() as f:
            test_artifact = os.path.join(TEST_DATA_DIR, "cat_tool_url.cwl")
            test_command = [
                "test",
                "--no-container",
                test_artifact,
            ]
            self._check_exit_code(test_command, exit_code=0)
            assert_exists(os.path.join(f, "tool_test_output.json"))

    @skip_if_environ("PLANEMO_SKIP_CWLTOOL_TESTS")
    @skip_unless_module("toil")
    def test_toil_tool_test(self):
        """Test testing a CWL tool with cwltool."""
        with self._isolate() as f:
            test_artifact = os.path.join(TEST_DATA_DIR, "int_tool.cwl")
            test_command = [
                "test",
                "--engine", "toil",
                test_artifact,
            ]
            self._check_exit_code(test_command, exit_code=0)
            assert_exists(os.path.join(f, "tool_test_output.json"))

    @skip_if_environ("PLANEMO_SKIP_CWLTOOL_TESTS")
    def test_output_checks_cwltool(self):
        """Test output assertions with a CWL tool with cwltool."""
        self._output_checks([])

    @skip_if_environ("PLANEMO_SKIP_CWLTOOL_TESTS")
    @skip_if_environ("PLANEMO_SKIP_SLOW_TESTS")
    @skip_unless_module("toil")
    def test_output_checks_toil(self):
        """Test output assertions with a CWL tool with toil."""
        self._output_checks(["--engine", "toil"])

    def _output_checks(self, extra_args):
        with self._isolate() as f:
            test_artifact = os.path.join(TEST_DATA_DIR, "output_tests_tool.cwl")
            test_command = [
                "test",
                "--no-container",
            ]
            test_command.extend(extra_args)
            test_command.append(test_artifact)
            self._check_exit_code(test_command, exit_code=1)
            output_json_path = os.path.join(f, "tool_test_output.json")
            with open(output_json_path, "r") as f:
                output = json.load(f)
            assert "tests" in output
            tests = output["tests"]
            # check out tests/data/output_tests_tool_test.yml
            expected_statuses = [
                "success",
                "failure",
                "success",
                "success",
                "failure",
                "success",
                "failure",
                "success",
                "failure",
                "success",
                "failure",
            ]
            for i in range(len(expected_statuses)):
                test_i = tests[i]
                data = test_i["data"]
                expected_status = expected_statuses[i]
                assert data["status"] == expected_status

    def append_profile_argument_if_needed(self, command):
        # Hook into tests to allow leveraging postgres databases to prevent Galaxy locking errors
        # while running tests.
        profile_name = os.getenv("PLANEMO_TEST_WORKFLOW_RUN_PROFILE", None)

        if profile_name:
            command += ["--profile", profile_name]

            database_type = os.getenv("PLANEMO_TEST_WORKFLOW_RUN_PROFILE_DATABASE_TYPE", None)
            if database_type:
                command += ["--database_type", database_type]

        return command
