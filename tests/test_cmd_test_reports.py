import json
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

    def test_markdown(self):
        with self._isolate() as f:
            json_path = os.path.join(TEST_DATA_DIR, "issue381.json")
            results_path = os.path.join(f, "markdown_results")
            self._check_exit_code(["test_reports", "--test_output_markdown", results_path, json_path], exit_code=0)
            assert os.path.exists(results_path)

            # Run minimal version
            minimal_results_path = os.path.join(f, "minimal_markdown_results")
            self._check_exit_code(
                ["test_reports", "--test_output_markdown_minimal", minimal_results_path, json_path], exit_code=0
            )
            assert os.path.exists(minimal_results_path)
            # Make sure minimal markdown is compacted
            assert os.path.getsize(minimal_results_path) < os.path.getsize(results_path)

    def test_markdown_counts_skip_status(self):
        with self._isolate() as f:
            json_path = os.path.join(f, "skipped.json")
            markdown_path = os.path.join(f, "markdown_results")
            report = {
                "version": "0.1",
                "tests": [
                    {
                        "id": "skipped-test",
                        "has_data": True,
                        "data": {"status": "skip"},
                    }
                ],
            }
            with open(json_path, "w") as handle:
                json.dump(report, handle)

            self._check_exit_code(["test_reports", "--test_output_markdown", markdown_path, json_path], exit_code=0)

            with open(markdown_path) as handle:
                markdown = handle.read()

            assert "| Skipped    | 1 |" in markdown

    def test_allure_normalizes_legacy_skipped_status(self):
        with self._isolate() as f:
            json_path = os.path.join(f, "skipped.json")
            results_path = os.path.join(f, "allure_results")
            report = {
                "version": "0.1",
                "tests": [
                    {
                        "id": "skipped-test",
                        "has_data": True,
                        "data": {"status": "skipped"},
                    }
                ],
            }
            with open(json_path, "w") as handle:
                json.dump(report, handle)

            self._check_exit_code(["test_reports", "--test_output_allure", results_path, json_path], exit_code=0)

            assert os.path.exists(results_path)

    def test_invalid_report_fails_with_friendly_error(self):
        with self._isolate() as f:
            json_path = os.path.join(f, "invalid.json")
            with open(json_path, "w") as handle:
                json.dump({"version": "0.1", "tests": [{"id": "bad", "has_data": True, "data": None}]}, handle)

            result = self._check_exit_code(["test_reports", json_path], exit_code=1)

            assert "Invalid Planemo test report JSON" in result.output
