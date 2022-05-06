import os

from .test_utils import (
    CliTestCase,
    TEST_DATA_DIR,
)


class CmdMergeReportsTestCase(CliTestCase):
    def test_merge_reports(self):
        with self._isolate():
            json_path = os.path.join(TEST_DATA_DIR, "issue381.json")
            self._check_exit_code(["merge_test_reports", json_path, json_path, json_path, "out.json"], exit_code=0)
