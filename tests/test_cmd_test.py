"""Module contains :class:`CmdTestTestCase` - integration tests for the ``test`` command."""
import json
import os
import shutil
from tempfile import (
    NamedTemporaryFile,
    TemporaryDirectory,
)

from .test_utils import (
    assert_exists,
    CliTestCase,
    PROJECT_TEMPLATES_DIR,
    run_verbosely,
    skip_if_environ,
    skip_unless_module,
    TEST_DATA_DIR,
    TEST_TOOLS_DIR,
)

FETCH_DATA_DATA_MANAGER_TEST_PATH = "data_manager/data_manager_fetch_genome_dbkeys_all_fasta/data_manager/data_manager_fetch_genome_all_fasta_dbkeys.xml"
BOWTIE2_DATA_MANAGER_TEST_PATH = (
    "data_manager/data_manager_bowtie2_index_builder/data_manager/bowtie2_index_builder.xml"
)


class CmdTestTestCase(CliTestCase):
    """Integration tests for the ``test`` command."""

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    def test_startup_timeout(self):
        """Test --galaxy_startup_timeout."""
        with self._isolate():
            test_artifact = os.path.join(TEST_DATA_DIR, FETCH_DATA_DATA_MANAGER_TEST_PATH)
            test_command = self._test_command(
                "--galaxy_startup_timeout", "1", test_artifact, "--no_dependency_resolution"
            )
            self._check_exit_code(test_command, exit_code=1)

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    def test_tool_in_directory(self):
        """Test with (single) tool in directory."""
        with self._isolate(), TemporaryDirectory() as tempdir:
            test_artifact = os.path.join(TEST_DATA_DIR, "tools", "ok_test_assert_command.xml")
            shutil.copy(test_artifact, tempdir)
            test_command = self._test_command(tempdir, "--no_dependency_resolution")
            self._check_exit_code(test_command, exit_code=0)

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    def test_data_manager(self):
        """Test testing a data manager test."""
        with self._isolate(), NamedTemporaryFile(prefix="data_manager_test_json") as json_out:
            test_artifact = os.path.join(TEST_DATA_DIR, FETCH_DATA_DATA_MANAGER_TEST_PATH)
            test_command = self._test_command("--test_output_json", json_out.name)
            test_command = self.append_profile_argument_if_needed(test_command)
            test_command += [
                "--no_dependency_resolution",
                test_artifact,
            ]
            self._check_exit_code(test_command, exit_code=0)
            with open(json_out.name) as fh:
                assert json.load(fh)["summary"]["num_tests"] == 1

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    def test_data_manager_docker_mount(self):
        """Test testing a data manager that needs (ro) access to the test-data folder."""
        with self._isolate(), NamedTemporaryFile(prefix="data_manager_test_json") as json_out:
            test_artifact = os.path.join(TEST_DATA_DIR, BOWTIE2_DATA_MANAGER_TEST_PATH)
            test_command = self._test_command("--test_output_json", json_out.name)
            test_command = self.append_profile_argument_if_needed(test_command)
            # data manager script is symlinked out of directory, will only work with `--docker_extra_volume`
            # we'll also add a bunch more to test multi path handling
            extra_volume = os.path.join(TEST_DATA_DIR, "data_manager")
            test_command += [
                "--no_dependency_resolution",
                "--biocontainers",
                "--docker_extra_volume",
                extra_volume,
                "--docker_extra_volume",
                extra_volume,
                test_artifact,
            ]
            self._check_exit_code(test_command, exit_code=0)
            with open(json_out.name) as fh:
                assert json.load(fh)["summary"]["num_tests"] == 1

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    def test_workflow_test_simple_yaml(self):
        """Test testing a simple YAML workflow with Galaxy."""
        with self._isolate():
            random_lines = os.path.join(PROJECT_TEMPLATES_DIR, "demo", "randomlines.xml")
            cat = os.path.join(PROJECT_TEMPLATES_DIR, "demo", "cat.xml")
            test_artifact = os.path.join(TEST_DATA_DIR, "wf1.gxwf.yml")
            test_command = self._test_command()
            test_command = self.append_profile_argument_if_needed(test_command)
            test_command += [
                "--no_dependency_resolution",
                "--extra_tools",
                random_lines,
                "--extra_tools",
                cat,
                test_artifact,
            ]
            self._check_exit_code(test_command, exit_code=0)

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    def test_tool_test_timeout(self):
        """Test if --test_timeout parameter is working."""
        with self._isolate(), NamedTemporaryFile(prefix="timeout_test_json") as json_out:
            test_artifact = os.path.join(TEST_TOOLS_DIR, "timeout.xml")
            test_command = self._test_command("--test_output_json", json_out.name)
            test_command = self.append_profile_argument_if_needed(test_command)
            test_command += [
                "--test_timeout",
                "1",
                test_artifact,
            ]
            self._check_exit_code(test_command, exit_code=1)
            with open(json_out.name) as fh:
                tool_test_json = json.load(fh)
                assert tool_test_json["summary"]["num_tests"] == 1
                # check run time, for smaller 10 since the test will take a bit longer than 1s
                # the important bit is that it's not about 30s (since the test tool calls `sleep 30`)
                assert (
                    float(tool_test_json["tests"][0]["data"]["time_seconds"]) <= 10
                ), "Test needed more than 10 sec but should time out after 1"
                assert (
                    "Timed out after" in tool_test_json["tests"][0]["data"]["output_problems"][0]
                ), "Time out did not happen"

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    def test_workflow_test_simple_yaml_dockerized(self):
        """Test testing a simple YAML workflow with Galaxy in Docker."""
        with self._isolate():
            random_lines = os.path.join(PROJECT_TEMPLATES_DIR, "demo", "randomlines.xml")
            cat = os.path.join(PROJECT_TEMPLATES_DIR, "demo", "cat.xml")
            test_artifact = os.path.join(TEST_DATA_DIR, "wf1.gxwf.yml")
            test_command = self._test_command()
            test_command = self.append_profile_argument_if_needed(test_command)
            test_command += [
                "--no_dependency_resolution",
                "--engine",
                "docker_galaxy",
                "--extra_tools",
                random_lines,
                "--extra_tools",
                cat,
                test_artifact,
            ]
            self._check_exit_code(test_command, exit_code=0)

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    def test_workflow_test_simple_ga(self):
        """Test testing a simple GA workflow with Galaxy."""
        with self._isolate():
            cat = os.path.join(PROJECT_TEMPLATES_DIR, "demo", "cat.xml")
            test_artifact = os.path.join(TEST_DATA_DIR, "wf2.ga")
            test_command = self._test_command()
            test_command = self.append_profile_argument_if_needed(test_command)
            test_command += [
                "--no_dependency_resolution",
                "--extra_tools",
                cat,
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
        """Test testing a workflow that uses distro tools."""
        with self._isolate():
            test_artifact = os.path.join(TEST_DATA_DIR, "wf4-distro-tools.gxwf.yml")
            test_command = self._test_command()
            test_command = self.append_profile_argument_if_needed(test_command)
            test_command += [
                "--no_dependency_resolution",
                test_artifact,
            ]
            # try:
            self._check_exit_code(test_command, exit_code=0)
            # except Exception:
            #    with open(os.path.join(f, "tool_test_output.json"), "r") as o:
            #        print(o.read())
            #    raise

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    def test_workflow_test_composite(self):
        """Test testing a workflow that uses composite inputs."""
        with self._isolate():
            test_artifact = os.path.join(TEST_DATA_DIR, "wf6-composite-inputs.gxwf.yml")
            composite_input_imzml = os.path.join(TEST_DATA_DIR, "composite_input_imzml.xml")
            test_command = self._test_command()
            test_command = self.append_profile_argument_if_needed(test_command)
            test_command += [
                "--no_dependency_resolution",
                "--extra_tools",
                composite_input_imzml,
                test_artifact,
            ]
            self._check_exit_code(test_command, exit_code=0)

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    def test_workflow_test_collection_inputs(self):
        """Test testing a workflow with collection inputs Galaxy."""
        with self._isolate():
            test_artifact = os.path.join(TEST_DATA_DIR, "wf5-collection-input.gxwf.yml")
            cat_list = os.path.join(TEST_DATA_DIR, "cat_list.xml")
            test_command = self._test_command()
            test_command = self.append_profile_argument_if_needed(test_command)
            test_command += [
                "--no_dependency_resolution",
                "--extra_tools",
                cat_list,
                test_artifact,
            ]
            # try:
            self._check_exit_code(test_command, exit_code=0)
            # except Exception:
            #    with open(os.path.join(f, "tool_test_output.json"), "r") as o:
            #        print(o.read())
            #    raise

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    def test_workflow_test_repository_installation_gxformat2(self):
        """Test testing a workflow with collection inputs Galaxy."""
        with self._isolate():
            test_artifact = os.path.join(TEST_DATA_DIR, "wf13_tool_shed_repository_gxformat2.yml")
            test_command = self._test_command()
            test_command = self.append_profile_argument_if_needed(test_command)
            test_command += [
                "--biocontainers",
                test_artifact,
            ]
            self._check_exit_code(test_command, exit_code=0)

    @skip_if_environ("PLANEMO_SKIP_CWLTOOL_TESTS")
    def test_cwltool_tool_test(self):
        """Test testing a CWL tool with cwltool."""
        with self._isolate() as f:
            test_artifact = os.path.join(TEST_DATA_DIR, "int_tool.cwl")
            test_command = self._test_command(test_artifact)
            self._check_exit_code(test_command, exit_code=0)
            assert_exists(os.path.join(f, "tool_test_output.json"))

    @skip_if_environ("PLANEMO_SKIP_CWLTOOL_TESTS")
    def test_cwltool_tool_url_inputs_test(self):
        """Test testing a CWL tool with cwltool."""
        with self._isolate() as f:
            test_artifact = os.path.join(TEST_DATA_DIR, "cat_tool_url.cwl")
            test_command = self._test_command(
                "--no-container",
                test_artifact,
            )
            self._check_exit_code(test_command, exit_code=0)
            assert_exists(os.path.join(f, "tool_test_output.json"))

    @skip_if_environ("PLANEMO_SKIP_CWLTOOL_TESTS")
    @skip_unless_module("toil")
    def test_toil_tool_test(self):
        """Test testing a CWL tool with cwltool."""
        with self._isolate() as f:
            test_artifact = os.path.join(TEST_DATA_DIR, "int_tool.cwl")
            test_command = self._test_command(
                "--engine",
                "toil",
                test_artifact,
            )
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

    def _test_command(self, *args):
        test_cmd = ["--verbose"] if run_verbosely() else []
        test_cmd.append("test")
        test_cmd.extend(args)
        return test_cmd

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
            with open(output_json_path) as f:
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

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    def test_workflow_with_optional_input_output_not_provided(self):
        with self._isolate():
            test_artifact = os.path.join(TEST_DATA_DIR, "wf16_optional_input_output_label.ga")
            test_command = self._test_command()
            test_command = self.append_profile_argument_if_needed(test_command)
            test_command.append(test_artifact)
            self._check_exit_code(test_command, exit_code=0)
