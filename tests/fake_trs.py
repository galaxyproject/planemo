"""Test fakes for TRS functionality."""

from typing import (
    Any,
)


class FakeTrsImporter:
    """Test fake that records imports without calling Galaxy."""

    def __init__(self, return_workflow: dict[str, Any] | None = None):
        """Initialize fake importer.

        Args:
            return_workflow: Workflow dict to return from imports (default: {"id": "fake_wf_id"})
        """
        self.imported_urls: list[str] = []
        self._return = return_workflow or {"id": "fake_wf_id"}

    def import_from_trs(self, trs_url: str) -> dict[str, Any]:
        """Record the import and return fake workflow.

        Args:
            trs_url: TRS URL being imported

        Returns:
            Fake workflow dict
        """
        self.imported_urls.append(trs_url)
        return self._return
