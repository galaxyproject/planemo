import json
import os

from planemo.test.models import PlanemoTestReport
from .test_utils import (
    CliTestCase,
    TEST_DATA_DIR,
)


class CmdMergeReportsTestCase(CliTestCase):
    def test_merge_reports(self):
        with self._isolate():
            json_path = os.path.join(TEST_DATA_DIR, "issue381.json")
            self._check_exit_code(["merge_test_reports", json_path, json_path, json_path, "out.json"], exit_code=0)

            with open("out.json") as f:
                merged_report = PlanemoTestReport.model_validate(json.load(f))

            assert merged_report.version == "0.1"
            assert merged_report.summary is not None
            assert merged_report.summary.num_tests == 12
            assert merged_report.exit_code == 1

    def test_merge_reports_normalizes_legacy_skipped_status(self):
        with self._isolate():
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
            with open("skipped.json", "w") as f:
                json.dump(report, f)

            self._check_exit_code(["merge_test_reports", "skipped.json", "out.json"], exit_code=0)

            with open("out.json") as f:
                merged_report = PlanemoTestReport.model_validate(json.load(f))

            assert merged_report.tests[0].data.status == "skip"
            assert merged_report.summary.num_skips == 1
            assert merged_report.exit_code == 0

    def test_merge_reports_invalid_report_fails_with_friendly_error(self):
        with self._isolate():
            with open("invalid.json", "w") as f:
                json.dump({"version": "0.1", "tests": [{"id": "bad", "has_data": True, "data": None}]}, f)

            result = self._check_exit_code(["merge_test_reports", "invalid.json", "out.json"], exit_code=1)

            assert "Invalid Planemo test report JSON" in result.output
            assert "invalid.json" in result.output
