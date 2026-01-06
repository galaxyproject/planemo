"""Tests for TRS ID resolution functionality."""

from unittest.mock import (
    Mock,
    patch,
)

from planemo.galaxy.workflows import (
    import_workflow_from_trs,
    parse_trs_id,
    parse_trs_uri,
    TRS_WORKFLOWS_PREFIX,
)
from planemo.runnable import RunnableType
from planemo.runnable_resolve import for_runnable_identifier


class TestTRSIdParsing:
    """Test TRS ID parsing to full URLs."""

    @patch("planemo.galaxy.workflows.requests.get")
    def test_parse_trs_id_workflow_without_version(self, mock_get):
        """Test parsing a workflow TRS ID without specific version fetches latest."""
        # Mock the Dockstore API response for versions
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.json.return_value = [
            {"name": "v0.2.0", "id": "version1"},
            {"name": "v0.1.14", "id": "version2"},
        ]
        mock_get.return_value = mock_response

        trs_id = "workflow/github.com/iwc-workflows/parallel-accession-download/main"
        result = parse_trs_id(trs_id)

        assert result is not None
        # Should fetch the first version from the list
        assert (
            result["trs_url"]
            == "https://dockstore.org/api/ga4gh/trs/v2/tools/#workflow/github.com/iwc-workflows/parallel-accession-download/main/versions/v0.2.0"
        )
        # Verify API was called to fetch versions
        mock_get.assert_called_once()

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

    @patch("planemo.galaxy.workflows.requests.get")
    def test_parse_trs_id_version_fetch_failure(self, mock_get):
        """Test parsing when version fetch fails falls back to default version."""
        # Mock a failed API request
        mock_get.side_effect = Exception("API error")

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

    @patch("planemo.galaxy.workflows.parse_trs_uri")
    def test_import_workflow_from_trs(self, mock_parse):
        """Test importing a workflow from TRS."""
        # Mock parse_trs_uri to return a full TRS URL
        expected_trs_url = "https://dockstore.org/api/ga4gh/trs/v2/tools/#workflow/github.com/org/repo/main"
        mock_parse.return_value = {"trs_url": expected_trs_url}

        # Mock Galaxy instance
        mock_gi = Mock()
        mock_workflows = Mock()
        mock_gi.workflows = mock_workflows

        # Mock the _make_url and _post methods
        mock_workflows._make_url.return_value = "https://galaxy.example.com/api/workflows"
        mock_workflows._post.return_value = {"id": "test_workflow_id", "name": "Test Workflow"}

        # Call import_workflow_from_trs
        trs_uri = "trs://workflow/github.com/org/repo/main"
        result = import_workflow_from_trs(trs_uri, mock_gi)

        # Verify the result
        assert result is not None
        assert result["id"] == "test_workflow_id"

        # Verify _post was called with correct payload
        mock_workflows._post.assert_called_once()
        call_args = mock_workflows._post.call_args
        assert call_args[1]["payload"]["trs_url"] == expected_trs_url

    @patch("planemo.galaxy.workflows.parse_trs_uri")
    def test_import_workflow_from_trs_with_version(self, mock_parse):
        """Test importing a workflow from TRS with specific version."""
        expected_trs_url = (
            "https://dockstore.org/api/ga4gh/trs/v2/tools/#workflow/github.com/org/repo/main/versions/v0.1.14"
        )
        mock_parse.return_value = {"trs_url": expected_trs_url}

        mock_gi = Mock()
        mock_workflows = Mock()
        mock_gi.workflows = mock_workflows
        mock_workflows._make_url.return_value = "https://galaxy.example.com/api/workflows"
        mock_workflows._post.return_value = {"id": "test_workflow_id"}

        trs_uri = "trs://workflow/github.com/org/repo/main/v0.1.14"
        result = import_workflow_from_trs(trs_uri, mock_gi)

        assert result is not None
        call_args = mock_workflows._post.call_args
        assert call_args[1]["payload"]["trs_url"] == expected_trs_url

    @patch("planemo.galaxy.workflows.parse_trs_uri")
    def test_import_workflow_from_trs_invalid_uri(self, mock_parse):
        """Test importing from an invalid TRS URI raises ValueError."""
        mock_parse.return_value = None

        mock_gi = Mock()
        trs_uri = "invalid_uri"

        try:
            import_workflow_from_trs(trs_uri, mock_gi)
            assert False, "Should have raised ValueError"
        except ValueError as e:
            assert "Invalid TRS URI" in str(e)


class TestTRSIdIntegration:
    """Test TRS ID integration with for_runnable_identifier."""

    @patch("planemo.runnable_resolve.translate_alias")
    def test_for_runnable_identifier_with_trs_id(self, mock_translate_alias):
        """Test that for_runnable_identifier creates TRS URIs."""
        # Mock translate_alias to return the input unchanged
        mock_translate_alias.side_effect = lambda ctx, identifier, profile: identifier

        ctx = Mock()
        trs_id = "workflow/github.com/iwc-workflows/parallel-accession-download/main"
        runnable = for_runnable_identifier(ctx, trs_id, {})

        assert runnable is not None
        assert runnable.type == RunnableType.galaxy_workflow
        assert runnable.uri == f"{TRS_WORKFLOWS_PREFIX}{trs_id}"
        assert runnable.is_trs_workflow_uri is True
        assert runnable.is_remote_workflow_uri is True

    @patch("planemo.runnable_resolve.translate_alias")
    def test_for_runnable_identifier_with_hash_prefix(self, mock_translate_alias):
        """Test that for_runnable_identifier handles # prefix."""
        mock_translate_alias.side_effect = lambda ctx, identifier, profile: identifier

        ctx = Mock()
        trs_id = "#workflow/github.com/iwc-workflows/parallel-accession-download/main/v0.1.14"
        runnable = for_runnable_identifier(ctx, trs_id, {})

        assert runnable is not None
        assert runnable.type == RunnableType.galaxy_workflow
        assert runnable.uri == f"{TRS_WORKFLOWS_PREFIX}{trs_id}"
        assert runnable.is_trs_workflow_uri is True

    @patch("planemo.runnable_resolve.translate_alias")
    def test_for_runnable_identifier_with_tool_trs_id(self, mock_translate_alias):
        """Test that for_runnable_identifier handles tool TRS IDs."""
        mock_translate_alias.side_effect = lambda ctx, identifier, profile: identifier

        ctx = Mock()
        trs_id = "tool/github.com/galaxyproject/example/v1.0"
        runnable = for_runnable_identifier(ctx, trs_id, {})

        assert runnable is not None
        assert runnable.uri == f"{TRS_WORKFLOWS_PREFIX}{trs_id}"
        assert runnable.is_trs_workflow_uri is True

    @patch("planemo.runnable_resolve.translate_alias")
    def test_for_runnable_identifier_with_full_dockstore_url(self, mock_translate_alias):
        """Test that for_runnable_identifier handles full Dockstore URLs."""
        mock_translate_alias.side_effect = lambda ctx, identifier, profile: identifier

        ctx = Mock()
        full_url = "https://dockstore.org/api/ga4gh/trs/v2/tools/#workflow/github.com/iwc-workflows/parallel-accession-download/main/versions/v0.1.14"
        runnable = for_runnable_identifier(ctx, full_url, {})

        assert runnable is not None
        assert runnable.type == RunnableType.galaxy_workflow
        # Should extract the path after the base URL
        expected_uri = "trs://#workflow/github.com/iwc-workflows/parallel-accession-download/main/versions/v0.1.14"
        assert runnable.uri == expected_uri
        assert runnable.is_trs_workflow_uri is True
