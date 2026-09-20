"""Tests for the planemo ``list_workflows`` command."""

import json
from unittest.mock import patch

import planemo.galaxy.api
from .test_utils import CliTestCase

WORKFLOWS = [
    {
        "id": "abc123",
        "name": "Published WF",
        "url": "/api/workflows/abc123",
        "source_metadata": {"url": "https://github.com/org/repo/wf.ga"},
        "published": True,
    },
    {
        "id": "def456",
        "name": "Unpublished WF",
        "url": "/api/workflows/def456",
        "source_metadata": None,
        "published": False,
    },
]


class FakeWorkflowClient:
    def get_workflows(self):
        return WORKFLOWS


class FakeGalaxyInstance:
    workflows = FakeWorkflowClient()


class CmdListWorkflowsTestCase(CliTestCase):
    def _list_workflows(self, galaxy_url, *extra_args):
        command = [
            "list_workflows",
            "--galaxy_url",
            galaxy_url,
            "--galaxy_user_key",
            "test_key",
            *extra_args,
        ]
        with patch.object(planemo.galaxy.api, "gi", lambda *args, **kwds: FakeGalaxyInstance()):
            return self._check_exit_code(command)

    def test_table_lists_each_workflow(self):
        result = self._list_workflows("http://localhost:8080")
        assert "Published WF" in result.output
        assert "Unpublished WF" in result.output
        assert "2 workflows found." in result.output

    def test_urls_are_browsable_not_api_endpoints(self):
        result = self._list_workflows("http://localhost:8080")
        # The API URL renders JSON - users need something they can open.
        assert "/api/workflows/abc123" not in result.output
        assert "http://localhost:8080/published/workflow?id=abc123" in result.output
        assert "http://localhost:8080/workflows/edit?id=def456" in result.output

    def test_trailing_slash_in_galaxy_url_does_not_double_up(self):
        result = self._list_workflows("http://localhost:8080/")
        assert "//published" not in result.output
        assert "//workflows" not in result.output
        assert "http://localhost:8080/published/workflow?id=abc123" in result.output

    def test_published_state_is_reported(self):
        result = self._list_workflows("http://localhost:8080")
        header = next(line for line in result.output.splitlines() if "Workflow ID" in line)
        assert "Published" in header
        published_row = next(line for line in result.output.splitlines() if line.startswith("abc123"))
        unpublished_row = next(line for line in result.output.splitlines() if line.startswith("def456"))
        assert published_row.split()[-3] == "yes"
        assert unpublished_row.split()[-2] == "no"

    def test_repo_url_shown_when_present_and_blank_otherwise(self):
        result = self._list_workflows("http://localhost:8080")
        assert "https://github.com/org/repo/wf.ga" in result.output
        assert "N/A" not in result.output

    def test_raw_emits_json_without_display_placeholders(self):
        result = self._list_workflows("http://localhost:8080", "--raw")
        parsed = json.loads(result.output)
        assert parsed["abc123"]["repo_url"] == "https://github.com/org/repo/wf.ga"
        # A missing repo URL must be null rather than a human-readable placeholder.
        assert parsed["def456"]["repo_url"] is None
        assert parsed["abc123"]["published"] is True

    def test_raw_suppresses_informational_output(self):
        result = self._list_workflows("http://localhost:8080", "--raw")
        assert "Looking for workflows" not in result.output
        assert "workflows found" not in result.output

    def test_missing_galaxy_url_reports_usage_error(self):
        result = self._runner.invoke(self._cli.planemo, ["list_workflows"])
        assert result.exit_code == 2
        assert "--galaxy_url" in result.output
        assert not isinstance(result.exception, AttributeError)
