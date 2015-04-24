import glob

from .test_utils import CliTestCase
from .test_utils import TEST_TOOLS_DIR


class LintTestCase(CliTestCase):

    def test_ok_tools(self):
        ok_tools = glob.glob("%s/ok_*" % TEST_TOOLS_DIR)
        for ok_tool in ok_tools:
            lint_cmd = ["lint", ok_tool]
            self._check_exit_code(lint_cmd)

    def test_fail_tools(self):
        fail_tools = glob.glob("%s/fail_*" % TEST_TOOLS_DIR)
        for fail_tool in fail_tools:
            lint_cmd = ["lint", fail_tool]
            self._check_exit_code(lint_cmd, exit_code=1)
