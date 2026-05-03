#!/usr/bin/env python

import json
import os

from planemo.test.models import PlanemoTestReport

from .test_utils import (
    CliTestCase,
    TEST_DATA_DIR,
)


class TestOutputSchema(CliTestCase):
    def test_output_schema_exports_test_report_schema(self):
        result = self._check_exit_code(["output_schema", "--format", "json", "--schema", "test-report"])
        metadata = json.loads(result.output)

        assert metadata["schema_version"] == "0.1"
        assert list(metadata["schemas"]) == ["test-report"]
        assert metadata["schemas"]["test-report"]["$schema"] == "https://json-schema.org/draft/2020-12/schema"

    def test_output_schema_exports_run_outputs_schema(self):
        result = self._check_exit_code(["output_schema", "--format", "json", "--schema", "run-outputs"])
        metadata = json.loads(result.output)

        assert list(metadata["schemas"]) == ["run-outputs"]
        assert metadata["schemas"]["run-outputs"]["$schema"] == "https://json-schema.org/draft/2020-12/schema"

    def test_issue381_validates_as_test_report(self):
        path = os.path.join(TEST_DATA_DIR, "issue381.json")
        with open(path) as f:
            report = json.load(f)

        parsed_report = PlanemoTestReport.model_validate(report)

        assert parsed_report.summary is not None
        assert parsed_report.summary.num_tests == 4
        assert parsed_report.tests[2].has_data is False
        assert parsed_report.tests[2].data is None

    def test_legacy_skipped_status_normalizes_to_skip(self):
        parsed_report = PlanemoTestReport.model_validate(
            {
                "version": "0.1",
                "tests": [
                    {
                        "id": "example",
                        "has_data": True,
                        "data": {"status": "skipped"},
                    }
                ],
            }
        )

        assert parsed_report.tests[0].data.status == "skip"
