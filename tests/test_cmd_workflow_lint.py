"""Tests for the ``workflow_lint`` command."""
import glob
import os

from .test_utils import (
    CliTestCase,
    TEST_DATA_DIR,
)


class CmdWorkflowLintTestCase(CliTestCase):
    def test_gxformat2_examples_as_repos(self):
        repos = glob.glob(_wf_repo("from_format2") + "/*")
        for repo in repos:
            repo_basename = os.path.basename(repo)
            try:
                expected_exit_code = int(repo_basename[0])
            except ValueError:
                # not a repo, just skip.
                continue
            lint_cmd = ["workflow_lint", "--skip", "tests,best_practices", repo]
            self._check_exit_code(lint_cmd, exit_code=expected_exit_code)

    def test_fail_level(self):
        # ensure missing tests normally cause it to fail...
        repo = _wf_repo("from_format2/0_basic_format2")
        lint_cmd = ["workflow_lint", repo]
        self._check_exit_code(lint_cmd, exit_code=1)

        # ... but not if fail_level is error
        repo = _wf_repo("from_format2/0_basic_format2")
        lint_cmd = ["workflow_lint", "--fail_level", "error", repo]
        self._check_exit_code(lint_cmd, exit_code=0)

    def test_workflow_test_linting(self):
        repo = _wf_repo("basic_format2_ok")
        lint_cmd = ["workflow_lint", "--skip", "best_practices", repo]
        self._check_exit_code(lint_cmd, exit_code=0)

        repo = _wf_repo("basic_native_ok")
        lint_cmd = ["workflow_lint", "--skip", "best_practices", repo]
        self._check_exit_code(lint_cmd, exit_code=0)

        repo = _wf_repo("basic_format2_input_misspelled")
        lint_cmd = ["workflow_lint", "--skip", "best_practices", repo]
        self._check_exit_code(lint_cmd, exit_code=1)

        repo = _wf_repo("basic_format2_input_missing")
        lint_cmd = ["workflow_lint", "--skip", "best_practices", repo]
        self._check_exit_code(lint_cmd, exit_code=1)

        repo = _wf_repo("basic_format2_output_misnamed")
        lint_cmd = ["workflow_lint", "--skip", "best_practices", repo]
        self._check_exit_code(lint_cmd, exit_code=1)

        repo = _wf_repo("basic_format2_missing_input")
        lint_cmd = ["workflow_lint", "--skip", "best_practices", repo]
        self._check_exit_code(lint_cmd, exit_code=1)

    def test_workflow_test_linting_control(self):
        # if we skip the tests linting - the above failures should pass
        repo = _wf_repo("basic_format2_input_misspelled")
        lint_cmd = ["workflow_lint", "--skip", "tests,best_practices", repo]
        self._check_exit_code(lint_cmd, exit_code=0)

        repo = _wf_repo("basic_format2_input_missing")
        lint_cmd = ["workflow_lint", "--skip", "tests,best_practices", repo]
        self._check_exit_code(lint_cmd, exit_code=0)

        repo = _wf_repo("basic_format2_output_misnamed")
        lint_cmd = ["workflow_lint", "--skip", "tests,best_practices", repo]
        self._check_exit_code(lint_cmd, exit_code=0)

        repo = _wf_repo("basic_format2_missing_input")
        lint_cmd = ["workflow_lint", "--skip", "tests,best_practices", repo]
        self._check_exit_code(lint_cmd, exit_code=0)

    def test_workflow_dockstore_linting(self):
        repo = _wf_repo("basic_format2_dockstore")
        lint_cmd = ["workflow_lint", "--skip", "best_practices", repo]
        self._check_exit_code(lint_cmd, exit_code=0)

        repo = _wf_repo("basic_format2_dockstore_empty")
        lint_cmd = ["workflow_lint", "--skip", "best_practices", repo]
        self._check_exit_code(lint_cmd, exit_code=1)

        repo = _wf_repo("basic_format2_dockstore_invalid_yaml")
        lint_cmd = ["workflow_lint", "--skip", "best_practices", repo]
        self._check_exit_code(lint_cmd, exit_code=1)

        repo = _wf_repo("basic_format2_dockstore_wrong_descriptor")
        lint_cmd = ["workflow_lint", "--skip", "best_practices", repo]
        self._check_exit_code(lint_cmd, exit_code=1)

        repo = _wf_repo("basic_format2_dockstore_wrong_test_file")
        lint_cmd = ["workflow_lint", "--skip", "best_practices", repo]
        self._check_exit_code(lint_cmd, exit_code=1)

        repo = _wf_repo("basic_format2_dockstore_misspelled_primary_key")
        lint_cmd = ["workflow_lint", "--skip", "best_practices", repo]
        self._check_exit_code(lint_cmd, exit_code=1)

    def test_workflow_dockstore_linting_control(self):
        # run same tests as above but make sure if we skip dockstore they
        # all pass
        repo = _wf_repo("basic_format2_dockstore")
        lint_cmd = ["workflow_lint", "--skip", "dockstore,best_practices", repo]
        self._check_exit_code(lint_cmd, exit_code=0)

        repo = _wf_repo("basic_format2_dockstore_empty")
        lint_cmd = ["workflow_lint", "--skip", "dockstore,best_practices", repo]
        self._check_exit_code(lint_cmd, exit_code=0)

        repo = _wf_repo("basic_format2_dockstore_invalid_yaml")
        lint_cmd = ["workflow_lint", "--skip", "dockstore,best_practices", repo]
        self._check_exit_code(lint_cmd, exit_code=0)

        repo = _wf_repo("basic_format2_dockstore_wrong_descriptor")
        lint_cmd = ["workflow_lint", "--skip", "dockstore,best_practices", repo]
        self._check_exit_code(lint_cmd, exit_code=0)

        repo = _wf_repo("basic_format2_dockstore_wrong_test_file")
        lint_cmd = ["workflow_lint", "--skip", "dockstore,best_practices", repo]
        self._check_exit_code(lint_cmd, exit_code=0)

        repo = _wf_repo("basic_format2_dockstore_misspelled_primary_key")
        lint_cmd = ["workflow_lint", "--skip", "dockstore,best_practices", repo]
        self._check_exit_code(lint_cmd, exit_code=0)

    def test_lint_test_examples(self):
        tags_wf = os.path.join(TEST_DATA_DIR, "wf10-tags-and-rules.gxwf.yml")
        lint_cmd = ["workflow_lint", "--skip", "best_practices", tags_wf]
        self._check_exit_code(lint_cmd, exit_code=0)

    def test_best_practices_linting_gx(self):
        workflow_path = "/".join((TEST_DATA_DIR, "wf14-unlinted-best-practices.yml"))
        lint_cmd = ["workflow_lint", workflow_path]
        result = self._runner.invoke(self._cli.planemo, lint_cmd)

        warnings = [
            "Workflow is not annotated.",
            "Workflow does not specify a creator.",
            "Workflow does not specify a license.",
            "Workflow step with ID None has no annotation.",
            "Workflow step with ID None has no label.",
            "Workflow missing test cases.",
            "Workflow step with ID None specifies an untyped parameter as an input.",
            "Workflow step with ID None specifies an untyped parameter in the post-job actions.",
        ]

        for warning in warnings:
            assert warning in result.output

    def test_best_practices_linting_ga(self):
        workflow_path = "/".join((TEST_DATA_DIR, "wf14-unlinted-best-practices.ga"))
        lint_cmd = ["workflow_lint", workflow_path]
        result = self._runner.invoke(self._cli.planemo, lint_cmd)

        warnings = [
            "Workflow is not annotated.",
            "Workflow does not specify a creator.",
            "Workflow does not specify a license.",
            "Workflow step with ID 0 has no annotation.",
            "Workflow step with ID 0 has no label.",
            "Workflow missing test cases.",
            "Workflow step with ID 1 specifies an untyped parameter as an input.",
            "Workflow step with ID 1 specifies an untyped parameter in the post-job actions.",
        ]

        for warning in warnings:
            assert warning in result.output

    def test_assertion_linting(self):
        workflow_path = "/".join((TEST_DATA_DIR, "wf15-test-assertions.yml"))
        lint_cmd = ["workflow_lint", workflow_path]
        result = self._runner.invoke(self._cli.planemo, lint_cmd)
        assert (
            "Invalid assertion in tests: assert_has_text got an unexpected keyword argument 'non_existent_attribute'"
            in result.output
        )

    def test_tool_id_linting_wrong_version(self):
        workflow_path = "/".join(
            (TEST_DATA_DIR, "wf_repos", "autoupdate_tests", "workflow_with_unexisting_version_of_tool.ga")
        )
        lint_cmd = ["workflow_lint", workflow_path]
        result = self._runner.invoke(self._cli.planemo, lint_cmd)
        assert (
            "ERROR: The tool toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_head_tool/0.1.0 is not in the toolshed"
            in result.output
        )

    def test_tool_id_linting_wrong_tool(self):
        workflow_path = "/".join((TEST_DATA_DIR, "wf_repos", "autoupdate_tests", "workflow_with_unexisting_tool.ga"))
        lint_cmd = ["workflow_lint", workflow_path]
        result = self._runner.invoke(self._cli.planemo, lint_cmd)
        assert "ERROR: The ToolShed returned an error when searching" in result.output

    def test_workflow_linting_asserts(self):
        repo = _wf_repo("basic_format2_ok_collection")
        lint_cmd = ["workflow_lint", "--skip", "best_practices", repo]
        self._check_exit_code(lint_cmd, exit_code=0)

        repo = _wf_repo("basic_format2_ok_list")
        lint_cmd = ["workflow_lint", "--skip", "best_practices", repo]
        self._check_exit_code(lint_cmd, exit_code=0)

        repo = _wf_repo("basic_format2_wrong_assert_list")
        lint_cmd = ["workflow_lint", "--skip", "best_practices", repo]
        result = self._runner.invoke(self._cli.planemo, lint_cmd)
        assert "ERROR: Invalid assertion in tests: assert_has_text missing a required argument: 'text'" in result.output

        repo = _wf_repo("basic_format2_collection_wrong_assert_list")
        lint_cmd = ["workflow_lint", "--skip", "best_practices", repo]
        result = self._runner.invoke(self._cli.planemo, lint_cmd)
        assert "ERROR: Invalid assertion in tests: assert_has_line missing a required argument: 'line'" in result.output

        repo = _wf_repo("basic_format2_collection_wrong_assert")
        lint_cmd = ["workflow_lint", "--skip", "best_practices", repo]
        result = self._runner.invoke(self._cli.planemo, lint_cmd)
        assert "ERROR: Invalid assertion in tests: assert_has_line missing a required argument: 'line'" in result.output


def _wf_repo(rel_path):
    return os.path.join(TEST_DATA_DIR, "wf_repos", rel_path)
