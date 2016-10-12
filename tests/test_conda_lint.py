"""Tests for the ``conda_lint`` command."""
import glob

from .test_utils import (
    CliTestCase,
    TEST_RECIPES_DIR,
)


class CondaLintTestCase(CliTestCase):
    """Container class defining test cases for the ``conda_lint`` command."""

    def test_ok_tools(self):
        """Iterate through the test recipes and check that ok tools are good."""
        ok_tools = glob.glob("%s/ok_*" % TEST_RECIPES_DIR)
        for ok_tool in ok_tools:
            lint_cmd = ["conda_lint", ok_tool]
            self._check_exit_code(lint_cmd)

    def test_fail_tools(self):
        """Iterate through the test recipes and check that failing tools fail."""
        fail_tools = glob.glob("%s/fail_*" % TEST_RECIPES_DIR)
        for fail_tool in fail_tools:
            lint_cmd = ["conda_lint", fail_tool]
            self._check_exit_code(lint_cmd, exit_code=1)
