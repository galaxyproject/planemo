import os

from .test_utils import CliTestCase, TEST_DATA_DIR


class CmdTestReportsTestCase(CliTestCase):

    def test_build_reports(self):
        with self._isolate():
            json_path = os.path.join(TEST_DATA_DIR, "issue381.json")
            self._check_exit_code(["test_reports", json_path], exit_code=0)
