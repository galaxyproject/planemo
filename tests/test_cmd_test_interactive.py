"""Tests for the ``test_interactive`` command."""

import json
import os
from unittest import mock

from .test_utils import (
    CliTestCase,
    mark,
    NON_ZERO_EXIT_CODE,
    skip_if_environ,
    skip_unless_executable,
    TEST_TOOLS_DIR,
)

INTERACTIVE_TOOL = os.path.join(TEST_TOOLS_DIR, "interactivetool_simple.xml")
INTERACTIVE_TOOL_WITH_INPUTS = os.path.join(TEST_TOOLS_DIR, "interactivetool_with_inputs.xml")
REGULAR_TOOL = os.path.join(TEST_TOOLS_DIR, "bwa_without_requirements.xml")

TEST_RUNNABLES = "planemo.commands.cmd_test_interactive.test_interactive_runnables"


class CmdTestInteractiveTestCase(CliTestCase):
    def test_no_interactive_tools_fails(self):
        # Only a regular tool is supplied - nothing to test.
        self._check_exit_code(["test_interactive", REGULAR_TOOL], exit_code=NON_ZERO_EXIT_CODE)

    def test_no_interactive_tools_message(self):
        result = self._check_exit_code(["test_interactive", REGULAR_TOOL], exit_code=NON_ZERO_EXIT_CODE)
        assert "No interactive tools found" in result.output

    def test_invalid_inputs_json_fails(self):
        self._check_exit_code(
            ["test_interactive", "--inputs", "{not valid json", INTERACTIVE_TOOL],
            exit_code=NON_ZERO_EXIT_CODE,
        )

    def test_inputs_must_be_object(self):
        self._check_exit_code(
            ["test_interactive", "--inputs", "[1, 2]", INTERACTIVE_TOOL], exit_code=NON_ZERO_EXIT_CODE
        )

    def test_unsupported_engine_fails(self):
        self._check_exit_code(
            ["test_interactive", "--engine", "cwltool", INTERACTIVE_TOOL], exit_code=NON_ZERO_EXIT_CODE
        )

    def test_happy_path_forwards_filtered_runnables_and_inputs(self):
        with self._isolate(), mock.patch(TEST_RUNNABLES) as mock_run:
            mock_run.return_value = 0
            self._check_exit_code(
                ["test_interactive", "--inputs", '{"param": "value"}', INTERACTIVE_TOOL, REGULAR_TOOL]
            )
        assert mock_run.call_count == 1
        # Positional args are (ctx, runnables); the runnables list is the second.
        runnables = mock_run.call_args.args[1]
        kwargs = mock_run.call_args.kwargs
        assert len(runnables) == 1
        assert runnables[0].path.endswith(os.path.join("tools", "interactivetool_simple.xml"))
        assert os.path.isfile(runnables[0].path)
        assert kwargs["inputs"] == {"param": "value"}
        assert kwargs["engine"] == "galaxy"
        assert kwargs["expect_content"] is None

    def test_expect_content_forwarded(self):
        with self._isolate(), mock.patch(TEST_RUNNABLES) as mock_run:
            mock_run.return_value = 0
            self._check_exit_code(["test_interactive", "--expect_content", "moo cow\n", INTERACTIVE_TOOL])
        kwargs = mock_run.call_args.kwargs
        assert kwargs["expect_content"] == "moo cow\n"

    def test_test_index_out_of_range_fails(self):
        # The tool has 2 test cases; --test_index 3 is invalid and must fail
        # fast, before any managed Galaxy is started.
        result = self._check_exit_code(
            ["test_interactive", "--test_index", "3", INTERACTIVE_TOOL_WITH_INPUTS],
            exit_code=NON_ZERO_EXIT_CODE,
        )
        assert "out of range" in result.output

    def test_test_index_no_tests_fails(self):
        # The simple tool declares no <tests>; --test_index is invalid for it.
        result = self._check_exit_code(
            ["test_interactive", "--test_index", "1", INTERACTIVE_TOOL],
            exit_code=NON_ZERO_EXIT_CODE,
        )
        assert "0 test case" in result.output

    def test_test_index_forwarded(self):
        with self._isolate(), mock.patch(TEST_RUNNABLES) as mock_run:
            mock_run.return_value = 0
            self._check_exit_code(["test_interactive", "--test_index", "1", INTERACTIVE_TOOL_WITH_INPUTS])
        assert mock_run.call_count == 1
        kwargs = mock_run.call_args.kwargs
        assert kwargs["test_index"] == 1
        runnables = mock_run.call_args.args[1]
        assert runnables[0].path.endswith(os.path.join("tools", "interactivetool_with_inputs.xml"))

    def test_test_index_and_inputs_both_forwarded(self):
        with self._isolate(), mock.patch(TEST_RUNNABLES) as mock_run:
            mock_run.return_value = 0
            self._check_exit_code(
                ["test_interactive", "--test_index", "1", "--inputs", '{"extra": "1"}', INTERACTIVE_TOOL_WITH_INPUTS]
            )
        kwargs = mock_run.call_args.kwargs
        assert kwargs["test_index"] == 1
        assert kwargs["inputs"] == {"extra": "1"}

    def test_exit_code_propagated(self):
        with self._isolate(), mock.patch(TEST_RUNNABLES) as mock_run:
            mock_run.return_value = 1
            self._check_exit_code(["test_interactive", INTERACTIVE_TOOL], exit_code=1)

    def test_docker_engine_allowed(self):
        with self._isolate(), mock.patch(TEST_RUNNABLES) as mock_run:
            mock_run.return_value = 0
            self._check_exit_code(["test_interactive", "--engine", "docker_galaxy", INTERACTIVE_TOOL])
        kwargs = mock_run.call_args.kwargs
        assert kwargs["engine"] == "docker_galaxy"


class TestInteractiveEngineResolvesInputs(CliTestCase):
    """The engine must resolve per-tool inputs (from --test_index / --inputs)
    and pass them to the tool run."""

    def _run_engine(self, **kwds):
        from planemo.engine.galaxy import LocalManagedGalaxyEngine
        from planemo.runnable import for_path

        runnable = for_path(INTERACTIVE_TOOL_WITH_INPUTS)
        engine = LocalManagedGalaxyEngine(ctx=mock.Mock())
        fake_config = mock.Mock()
        success_case = {
            "id": "interactivetool_with_inputs-interactive",
            "has_data": True,
            "test_type": "interactive",
            "data": {"status": "success"},
        }
        with (
            mock.patch.object(LocalManagedGalaxyEngine, "ensure_runnables_served") as mock_served,
            mock.patch("planemo.engine.galaxy.run_interactive_tool_test") as mock_run,
        ):
            mock_served.return_value.__enter__.return_value = fake_config
            mock_run.return_value = success_case
            engine.test_interactive([runnable], **kwds)
        return mock_run

    def test_test_index_inputs_reach_tool_run(self):
        mock_run = self._run_engine(test_index=1)
        assert mock_run.call_count == 1
        # run_interactive_tool_test(config, tool_id, inputs=...)
        assert mock_run.call_args.args[1] == "interactivetool_with_inputs"
        assert mock_run.call_args.kwargs["inputs"] == {"greeting": "moo cow"}

    def test_manual_inputs_override_test_index(self):
        mock_run = self._run_engine(test_index=1, inputs={"greeting": "override"})
        assert mock_run.call_args.kwargs["inputs"] == {"greeting": "override"}

    def test_no_inputs_yields_empty(self):
        mock_run = self._run_engine()
        assert mock_run.call_args.kwargs["inputs"] == {}


class CmdTestInteractiveIntegrationTestCase(CliTestCase):
    """End-to-end test that runs ``test_interactive`` against a managed Galaxy.

    This launches a real (disposable) Galaxy instance and runs the sample
    interactive tool in a Docker container. It is gated the same way as
    :meth:`ServeTestCase.test_serve_interactivetool` and requires a container
    runtime to be available.
    """

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    @skip_unless_executable("docker")
    @mark.tests_galaxy_branch
    def test_interactive_tool_startup_and_content(self):
        with self._isolate():
            result = self._check_exit_code(
                [
                    "test_interactive",
                    "--biocontainers",
                    "--expect_content",
                    "moo cow\n",
                    INTERACTIVE_TOOL,
                ]
            )
            assert "passed" in result.output
            with open("tool_test_output.json") as f:
                report = json.load(f)
            summary = report["summary"]
            assert summary["num_tests"] == 1
            assert summary["num_failures"] == 0
            assert summary["num_errors"] == 0
            test_case = report["tests"][0]
            assert test_case["test_type"] == "interactive"
            assert test_case["data"]["status"] == "success"
            assert test_case["data"]["content"] == "moo cow\n"

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    @skip_unless_executable("docker")
    @mark.tests_galaxy_branch
    def test_interactive_tool_content_mismatch_fails(self):
        with self._isolate():
            result = self._check_exit_code(
                [
                    "test_interactive",
                    "--biocontainers",
                    "--expect_content",
                    "not the right content\n",
                    INTERACTIVE_TOOL,
                ],
                exit_code=NON_ZERO_EXIT_CODE,
            )
            assert "failed" in result.output
            with open("tool_test_output.json") as f:
                report = json.load(f)
            test_case = report["tests"][0]
            assert test_case["data"]["status"] == "failure"
            assert any("does not match" in problem for problem in test_case["data"]["output_problems"])

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    @skip_unless_executable("docker")
    @mark.tests_galaxy_branch
    def test_interactive_tool_requiring_params_runs_test_case(self):
        # An interactive tool that declares a required parameter cannot start
        # unless it is given one. Running it with --test_index supplies the
        # test case's <param> as the tool input so the tool actually launches.
        with self._isolate():
            result = self._check_exit_code(
                [
                    "test_interactive",
                    "--biocontainers",
                    "--test_index",
                    "1",
                    "--expect_content",
                    "moo cow\n",
                    INTERACTIVE_TOOL_WITH_INPUTS,
                ]
            )
            assert "passed" in result.output
            with open("tool_test_output.json") as f:
                report = json.load(f)
            test_case = report["tests"][0]
            assert test_case["data"]["status"] == "success"
            # The inputs came from test case 1 of the tool's <tests>.
            assert test_case["data"]["inputs"] == {"greeting": "moo cow"}
            assert test_case["data"]["content"] == "moo cow\n"
