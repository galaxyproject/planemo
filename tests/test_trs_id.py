"""Tests for TRS ID resolution functionality."""

import json
from pathlib import Path
from unittest.mock import (
    Mock,
    patch,
)

import responses

import planemo.cli
from planemo.galaxy.workflows import (
    import_workflow_from_trs,
    parse_trs_id,
    parse_trs_uri,
    TRS_WORKFLOWS_PREFIX,
)
from planemo.runnable import RunnableType
from planemo.runnable_resolve import for_runnable_identifier
from .fake_trs import FakeTrsImporter

FIXTURES = Path(__file__).parent / "fixtures" / "dockstore"


class TestTRSIdParsing:
    """Test TRS ID parsing to full URLs."""

    @responses.activate
    def test_parse_trs_id_workflow_without_version(self):
        """Test parsing a workflow TRS ID without specific version fetches latest."""
        # Load fixture data from recorded API response
        fixture = json.loads((FIXTURES / "versions_iwc_parallel_accession.json").read_text())

        responses.add(
            responses.GET,
            "https://dockstore.org/api/ga4gh/trs/v2/tools/%23workflow%2Fgithub.com%2Fiwc-workflows%2Fparallel-accession-download%2Fmain/versions",
            json=fixture,
            status=200,
        )

        trs_id = "workflow/github.com/iwc-workflows/parallel-accession-download/main"
        result = parse_trs_id(trs_id)

        assert result is not None
        # Should fetch the first version from the list
        first_version = fixture[0]["name"]
        assert (
            result["trs_url"]
            == f"https://dockstore.org/api/ga4gh/trs/v2/tools/#workflow/github.com/iwc-workflows/parallel-accession-download/main/versions/{first_version}"
        )

    def test_parse_trs_id_workflow_with_version(self):
        """Test parsing a workflow TRS ID with specific version."""
        trs_id = "workflow/github.com/iwc-workflows/parallel-accession-download/main/v0.1.14"
        result = parse_trs_id(trs_id)

        assert result is not None
        assert (
            result["trs_url"]
            == "https://dockstore.org/api/ga4gh/trs/v2/tools/#workflow/github.com/iwc-workflows/parallel-accession-download/main/versions/v0.1.14"
        )

    def test_parse_trs_id_with_hash_prefix(self):
        """Test parsing a TRS ID with # prefix."""
        trs_id = "#workflow/github.com/iwc-workflows/parallel-accession-download/main/v0.1.14"
        result = parse_trs_id(trs_id)

        assert result is not None
        assert (
            result["trs_url"]
            == "https://dockstore.org/api/ga4gh/trs/v2/tools/#workflow/github.com/iwc-workflows/parallel-accession-download/main/versions/v0.1.14"
        )

    def test_parse_trs_id_tool(self):
        """Test parsing a tool TRS ID."""
        trs_id = "tool/github.com/galaxyproject/example-tool/main/v1.0"
        result = parse_trs_id(trs_id)

        assert result is not None
        assert (
            result["trs_url"]
            == "https://dockstore.org/api/ga4gh/trs/v2/tools/#tool/github.com/galaxyproject/example-tool/main/versions/v1.0"
        )

    @responses.activate
    def test_parse_trs_id_version_fetch_failure(self):
        """Test parsing when version fetch fails falls back to default version."""
        # Simulate a failed API request
        responses.add(
            responses.GET,
            "https://dockstore.org/api/ga4gh/trs/v2/tools/%23workflow%2Fgithub.com%2Forg%2Frepo%2Fmain/versions",
            json={"error": "Not found"},
            status=404,
        )

        trs_id = "workflow/github.com/org/repo/main"
        result = parse_trs_id(trs_id)

        assert result is not None
        # Should fallback to using workflow_name as default version (Galaxy requires /versions/)
        assert (
            result["trs_url"]
            == "https://dockstore.org/api/ga4gh/trs/v2/tools/#workflow/github.com/org/repo/main/versions/main"
        )

    def test_parse_trs_id_invalid(self):
        """Test parsing invalid TRS IDs."""
        assert parse_trs_id("invalid") is None
        assert parse_trs_id("workflow/github.com/org") is None  # Too few parts


class TestTRSUriParsing:
    """Test TRS URI parsing."""

    @patch("planemo.galaxy.workflows.requests.get")
    def test_parse_trs_uri_workflow(self, mock_get):
        """Test parsing a workflow TRS URI."""
        # Mock the Dockstore API response
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.json.return_value = [{"name": "v0.1.14"}]
        mock_get.return_value = mock_response

        trs_uri = "trs://workflow/github.com/iwc-workflows/parallel-accession-download/main"
        result = parse_trs_uri(trs_uri)

        assert result is not None
        assert (
            result["trs_url"]
            == "https://dockstore.org/api/ga4gh/trs/v2/tools/#workflow/github.com/iwc-workflows/parallel-accession-download/main/versions/v0.1.14"
        )

    def test_parse_trs_uri_with_version(self):
        """Test parsing a TRS URI with version."""
        trs_uri = "trs://workflow/github.com/org/repo/main/v0.1.14"
        result = parse_trs_uri(trs_uri)

        assert result is not None
        assert (
            result["trs_url"]
            == "https://dockstore.org/api/ga4gh/trs/v2/tools/#workflow/github.com/org/repo/main/versions/v0.1.14"
        )

    def test_parse_trs_uri_from_full_url(self):
        """Test parsing a TRS URI created from a full Dockstore URL."""
        # This simulates what happens when user provides a full URL
        trs_uri = "trs://#workflow/github.com/iwc-workflows/parallel-accession-download/main/versions/v0.1.14"
        result = parse_trs_uri(trs_uri)

        assert result is not None
        expected_url = "https://dockstore.org/api/ga4gh/trs/v2/tools/#workflow/github.com/iwc-workflows/parallel-accession-download/main/versions/v0.1.14"
        assert result["trs_url"] == expected_url

    def test_parse_trs_uri_invalid(self):
        """Test parsing an invalid TRS URI."""
        assert parse_trs_uri("invalid") is None
        assert parse_trs_uri("trs://invalid") is None
        assert parse_trs_uri("gxid://workflows/abc") is None


class TestTRSWorkflowImport:
    """Test TRS workflow import."""

    def test_import_workflow_from_trs(self):
        """Test importing a workflow from TRS."""
        fake = FakeTrsImporter(return_workflow={"id": "test_wf_id", "name": "Test Workflow"})
        trs_uri = "trs://workflow/github.com/org/repo/main/v0.1.14"

        result = import_workflow_from_trs(trs_uri, fake)

        assert result is not None
        assert result["id"] == "test_wf_id"
        assert result["name"] == "Test Workflow"
        # Verify the importer was called with the correct URL
        assert len(fake.imported_urls) == 1
        assert (
            fake.imported_urls[0]
            == "https://dockstore.org/api/ga4gh/trs/v2/tools/#workflow/github.com/org/repo/main/versions/v0.1.14"
        )

    def test_import_workflow_from_trs_with_version(self):
        """Test importing a workflow from TRS with specific version."""
        fake = FakeTrsImporter(return_workflow={"id": "test_wf_id"})
        trs_uri = "trs://workflow/github.com/org/repo/main/v0.1.14"

        result = import_workflow_from_trs(trs_uri, fake)

        assert result is not None
        assert len(fake.imported_urls) == 1
        assert (
            fake.imported_urls[0]
            == "https://dockstore.org/api/ga4gh/trs/v2/tools/#workflow/github.com/org/repo/main/versions/v0.1.14"
        )

    def test_import_workflow_from_trs_invalid_uri(self):
        """Test importing from an invalid TRS URI raises ValueError."""
        fake = FakeTrsImporter()
        trs_uri = "invalid_uri"

        try:
            import_workflow_from_trs(trs_uri, fake)
            assert False, "Should have raised ValueError"
        except ValueError as e:
            assert "Invalid TRS URI" in str(e)
        # Verify no imports were attempted
        assert len(fake.imported_urls) == 0


class TestTRSIdIntegration:
    """Test TRS ID integration with for_runnable_identifier."""

    def test_for_runnable_identifier_with_trs_id(self, tmp_path):
        """Test that for_runnable_identifier creates TRS URIs."""
        ctx = planemo.cli.PlanemoCliContext()
        ctx.planemo_directory = str(tmp_path)

        trs_id = "workflow/github.com/iwc-workflows/parallel-accession-download/main"
        runnable = for_runnable_identifier(ctx, trs_id, {})

        assert runnable is not None
        assert runnable.type == RunnableType.galaxy_workflow
        assert runnable.uri == f"{TRS_WORKFLOWS_PREFIX}{trs_id}"
        assert runnable.is_trs_workflow_uri is True
        assert runnable.is_remote_workflow_uri is True

    def test_for_runnable_identifier_with_hash_prefix(self, tmp_path):
        """Test that for_runnable_identifier handles # prefix."""
        ctx = planemo.cli.PlanemoCliContext()
        ctx.planemo_directory = str(tmp_path)

        trs_id = "#workflow/github.com/iwc-workflows/parallel-accession-download/main/v0.1.14"
        runnable = for_runnable_identifier(ctx, trs_id, {})

        assert runnable is not None
        assert runnable.type == RunnableType.galaxy_workflow
        assert runnable.uri == f"{TRS_WORKFLOWS_PREFIX}{trs_id}"
        assert runnable.is_trs_workflow_uri is True

    def test_for_runnable_identifier_with_tool_trs_id(self, tmp_path):
        """Test that for_runnable_identifier handles tool TRS IDs."""
        ctx = planemo.cli.PlanemoCliContext()
        ctx.planemo_directory = str(tmp_path)

        trs_id = "tool/github.com/galaxyproject/example/v1.0"
        runnable = for_runnable_identifier(ctx, trs_id, {})

        assert runnable is not None
        assert runnable.uri == f"{TRS_WORKFLOWS_PREFIX}{trs_id}"
        assert runnable.is_trs_workflow_uri is True

    def test_for_runnable_identifier_with_full_dockstore_url(self, tmp_path):
        """Test that for_runnable_identifier handles full Dockstore URLs."""
        ctx = planemo.cli.PlanemoCliContext()
        ctx.planemo_directory = str(tmp_path)

        full_url = "https://dockstore.org/api/ga4gh/trs/v2/tools/#workflow/github.com/iwc-workflows/parallel-accession-download/main/versions/v0.1.14"
        runnable = for_runnable_identifier(ctx, full_url, {})

        assert runnable is not None
        assert runnable.type == RunnableType.galaxy_workflow
        # Should extract the path after the base URL
        expected_uri = "trs://#workflow/github.com/iwc-workflows/parallel-accession-download/main/versions/v0.1.14"
        assert runnable.uri == expected_uri
        assert runnable.is_trs_workflow_uri is True
