import os

from .test_utils import (
    CliTestCase,
    TEST_DATA_DIR,
)


class CmdTestReportsTestCase(CliTestCase):
    def test_build_reports(self):
        with self._isolate():
            json_path = os.path.join(TEST_DATA_DIR, "issue381.json")
            self._check_exit_code(["test_reports", json_path], exit_code=0)

    def test_allure(self):
        with self._isolate() as f:
            json_path = os.path.join(TEST_DATA_DIR, "issue381.json")
            results_path = os.path.join(f, "allure_results")
            self._check_exit_code(["test_reports", "--test_output_allure", results_path, json_path], exit_code=0)
            assert os.path.exists(results_path)
            assert os.path.isdir(results_path)
            assert len(os.listdir(results_path))
