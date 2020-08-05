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
            lint_cmd = ["workflow_lint", "--skip", "tests", repo]
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
        lint_cmd = ["workflow_lint", repo]
        self._check_exit_code(lint_cmd, exit_code=0)

        repo = _wf_repo("basic_native_ok")
        lint_cmd = ["workflow_lint", repo]
        self._check_exit_code(lint_cmd, exit_code=0)

        repo = _wf_repo("basic_format2_input_misspelled")
        lint_cmd = ["workflow_lint", repo]
        self._check_exit_code(lint_cmd, exit_code=1)

        repo = _wf_repo("basic_format2_input_missing")
        lint_cmd = ["workflow_lint", repo]
        self._check_exit_code(lint_cmd, exit_code=1)

        repo = _wf_repo("basic_format2_output_misnamed")
        lint_cmd = ["workflow_lint", repo]
        self._check_exit_code(lint_cmd, exit_code=1)

        repo = _wf_repo("basic_format2_missing_input")
        lint_cmd = ["workflow_lint", repo]
        self._check_exit_code(lint_cmd, exit_code=1)

    def test_workflow_test_linting_control(self):
        # if we skip the tests linting - the above failures should pass
        repo = _wf_repo("basic_format2_input_misspelled")
        lint_cmd = ["workflow_lint", "--skip", "tests", repo]
        self._check_exit_code(lint_cmd, exit_code=0)

        repo = _wf_repo("basic_format2_input_missing")
        lint_cmd = ["workflow_lint", "--skip", "tests", repo]
        self._check_exit_code(lint_cmd, exit_code=0)

        repo = _wf_repo("basic_format2_output_misnamed")
        lint_cmd = ["workflow_lint", "--skip", "tests", repo]
        self._check_exit_code(lint_cmd, exit_code=0)

        repo = _wf_repo("basic_format2_missing_input")
        lint_cmd = ["workflow_lint", "--skip", "tests", repo]
        self._check_exit_code(lint_cmd, exit_code=0)

    def test_workflow_dockstore_linting(self):
        repo = _wf_repo("basic_format2_dockstore")
        lint_cmd = ["workflow_lint", repo]
        self._check_exit_code(lint_cmd, exit_code=0)

        repo = _wf_repo("basic_format2_dockstore_empty")
        lint_cmd = ["workflow_lint", repo]
        self._check_exit_code(lint_cmd, exit_code=1)

        repo = _wf_repo("basic_format2_dockstore_invalid_yaml")
        lint_cmd = ["workflow_lint", repo]
        self._check_exit_code(lint_cmd, exit_code=1)

        repo = _wf_repo("basic_format2_dockstore_wrong_descriptor")
        lint_cmd = ["workflow_lint", repo]
        self._check_exit_code(lint_cmd, exit_code=1)

        repo = _wf_repo("basic_format2_dockstore_wrong_test_file")
        lint_cmd = ["workflow_lint", repo]
        self._check_exit_code(lint_cmd, exit_code=1)

        repo = _wf_repo("basic_format2_dockstore_misspelled_primary_key")
        lint_cmd = ["workflow_lint", repo]
        self._check_exit_code(lint_cmd, exit_code=1)

    def test_workflow_dockstore_linting_control(self):
        # run same tests as above but make sure if we skip dockstore they
        # all pass
        repo = _wf_repo("basic_format2_dockstore")
        lint_cmd = ["workflow_lint", "--skip", "dockstore", repo]
        self._check_exit_code(lint_cmd, exit_code=0)

        repo = _wf_repo("basic_format2_dockstore_empty")
        lint_cmd = ["workflow_lint", "--skip", "dockstore", repo]
        self._check_exit_code(lint_cmd, exit_code=0)

        repo = _wf_repo("basic_format2_dockstore_invalid_yaml")
        lint_cmd = ["workflow_lint", "--skip", "dockstore", repo]
        self._check_exit_code(lint_cmd, exit_code=0)

        repo = _wf_repo("basic_format2_dockstore_wrong_descriptor")
        lint_cmd = ["workflow_lint", "--skip", "dockstore", repo]
        self._check_exit_code(lint_cmd, exit_code=0)

        repo = _wf_repo("basic_format2_dockstore_wrong_test_file")
        lint_cmd = ["workflow_lint", "--skip", "dockstore", repo]
        self._check_exit_code(lint_cmd, exit_code=0)

        repo = _wf_repo("basic_format2_dockstore_misspelled_primary_key")
        lint_cmd = ["workflow_lint", "--skip", "dockstore", repo]
        self._check_exit_code(lint_cmd, exit_code=0)

    def test_lint_test_examples(self):
        tags_wf = os.path.join(TEST_DATA_DIR, "wf10-tags-and-rules.gxwf.yml")
        lint_cmd = ["workflow_lint", tags_wf]
        self._check_exit_code(lint_cmd, exit_code=0)


def _wf_repo(rel_path):
    return os.path.join(TEST_DATA_DIR, "wf_repos", rel_path)
