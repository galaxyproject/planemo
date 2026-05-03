#!/usr/bin/env python

import json
import os

from planemo.commands.cmd_invocation_download import invocation_download_manifest
from planemo.output_models import PlanemoInvocationDownloadManifest
from planemo.galaxy.test.actions import handle_reports
from planemo.test.models import PlanemoTestReport

from .test_utils import (
    CliTestCase,
    TEST_DATA_DIR,
)


class DummyContext:
    def vlog(self, message, exception=False):
        pass


class TestOutputSchema(CliTestCase):
    def test_output_schema_exports_test_report_schema(self):
        result = self._check_exit_code(["output_schema", "--schema", "test-report"])
        metadata = json.loads(result.output)

        assert metadata["schema_version"] == "0.1"
        assert list(metadata["schemas"]) == ["test-report"]
        assert metadata["schemas"]["test-report"]["$schema"] == "https://json-schema.org/draft/2020-12/schema"

    def test_output_schema_exports_run_outputs_schema(self):
        result = self._check_exit_code(["output_schema", "--schema", "run-outputs"])
        metadata = json.loads(result.output)

        assert list(metadata["schemas"]) == ["run-outputs"]
        assert metadata["schemas"]["run-outputs"]["$schema"] == "https://json-schema.org/draft/2020-12/schema"

    def test_output_schema_exports_invocation_download_manifest_schema(self):
        result = self._check_exit_code(
            ["output_schema", "--schema", "invocation-download-manifest"]
        )
        metadata = json.loads(result.output)

        assert list(metadata["schemas"]) == ["invocation-download-manifest"]
        assert (
            metadata["schemas"]["invocation-download-manifest"]["$schema"]
            == "https://json-schema.org/draft/2020-12/schema"
        )

    def test_output_schema_exports_cli_metadata_schema(self):
        result = self._check_exit_code(["output_schema", "--schema", "cli-metadata"])
        metadata = json.loads(result.output)

        assert list(metadata["schemas"]) == ["cli-metadata"]
        assert metadata["schemas"]["cli-metadata"]["$schema"] == "https://json-schema.org/draft/2020-12/schema"

    def test_invocation_download_manifest_relative_paths(self):
        output_directory = os.path.join(os.getcwd(), "outputs")
        output_path = os.path.join(output_directory, "answer.txt")
        manifest = invocation_download_manifest(
            "invocation123",
            output_directory,
            {
                "answer": {"path": output_path, "basename": "answer.txt", "class": "File"},
                "optional": None,
            },
            output_json=os.path.join(output_directory, "manifest.json"),
            missing_output_reasons={"optional": "skipped"},
        )

        parsed_manifest = PlanemoInvocationDownloadManifest.model_validate(manifest)

        assert parsed_manifest.path_type == "relative"
        assert parsed_manifest.output_directory == "."
        assert parsed_manifest.outputs["answer"]["path"] == "answer.txt"
        assert parsed_manifest.output_json == "manifest.json"
        assert parsed_manifest.missing_outputs[0].id == "optional"
        assert parsed_manifest.missing_outputs[0].reason == "skipped"

    def test_invocation_download_manifest_absolute_paths(self):
        output_directory = os.path.join(os.getcwd(), "outputs")
        output_path = os.path.join(output_directory, "answer.txt")
        manifest = invocation_download_manifest(
            "invocation123",
            output_directory,
            {"answer": {"path": output_path, "basename": "answer.txt", "class": "File"}},
            output_json=os.path.join(output_directory, "manifest.json"),
            path_type="absolute",
        )

        parsed_manifest = PlanemoInvocationDownloadManifest.model_validate(manifest)

        assert parsed_manifest.path_type == "absolute"
        assert parsed_manifest.output_directory == output_directory
        assert parsed_manifest.outputs["answer"]["path"] == output_path
        assert parsed_manifest.output_json == os.path.join(output_directory, "manifest.json")

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

    def test_handle_reports_validates_test_output_json_before_writing(self):
        with self._isolate() as f:
            output_path = os.path.join(f, "invalid.json")
            invalid_report = {
                "version": "0.1",
                "tests": [
                    {
                        "id": "invalid",
                        "has_data": True,
                        "data": None,
                    }
                ],
            }

            try:
                handle_reports(DummyContext(), invalid_report, {"test_output_json": output_path})
            except ValueError:
                pass
            else:
                raise AssertionError("Invalid test report was written without validation error")

            assert not os.path.exists(output_path)

    def test_test_report_requires_status_when_data_present(self):
        try:
            PlanemoTestReport.model_validate(
                {
                    "version": "0.1",
                    "tests": [
                        {
                            "id": "invalid",
                            "has_data": True,
                            "data": {},
                        }
                    ],
                }
            )
        except ValueError:
            pass
        else:
            raise AssertionError("Test report with data but no status validated")
