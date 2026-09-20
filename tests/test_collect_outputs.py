"""Tests for output collection in :mod:`planemo.galaxy.activity`."""

from typing import (
    Any,
    cast,
    Dict,
    List,
    Optional,
    TYPE_CHECKING,
)

import pytest

from planemo.galaxy import activity
from planemo.galaxy.activity import GalaxyBaseRunResponse

if TYPE_CHECKING:
    from bioblend.galaxy import GalaxyInstance

    from planemo.cli import PlanemoCliContext

SAMPLE_SHEET_COLLECTION_METADATA: Dict[str, Any] = {
    "id": "e9d955ca403027d5",
    "model_class": "HistoryDatasetCollectionAssociation",
    "history_content_type": "dataset_collection",
    "collection_type": "sample_sheet:paired",
    "populated": True,
    "element_count": 1,
    "elements": [
        {
            "id": "dd871f1d42ebeefa",
            "model_class": "DatasetCollectionElement",
            "element_index": 0,
            "element_identifier": "sample1",
            "element_type": "dataset_collection",
            "object": {
                "id": "b1e4c6f0a2d3e5f7",
                "model_class": "DatasetCollection",
                "collection_type": "paired",
                "populated": True,
                "elements": [],
            },
        }
    ],
}


class MockContext:
    def log(self, msg, *args, **kwds):
        pass

    def vlog(self, msg, *args, **kwds):
        pass


class MockRunnableOutput:
    def __init__(self, output_id: str) -> None:
        self._id = output_id

    def get_id(self) -> str:
        return self._id

    def is_optional(self) -> bool:
        return False


class CollectionRunResponse(GalaxyBaseRunResponse):
    """A run response whose single output is an unsupported collection type."""

    def to_galaxy_output(self, runnable_output):
        return None

    def output_src(self, output, ignore_missing_outputs: Optional[bool] = False) -> Dict[str, str]:
        return {"src": "hdca", "id": "e9d955ca403027d5"}

    def _get_metadata(self, history_content_type: str, content_id: str) -> Dict[str, Any]:
        assert history_content_type == "dataset_collection"
        return dict(SAMPLE_SHEET_COLLECTION_METADATA)


@pytest.fixture
def unsupported_collection_response(monkeypatch) -> CollectionRunResponse:
    """A response where output_to_cwl_json returns None, as it does for sample_sheet."""

    def mock_get_outputs(runnable, gi=None) -> List[MockRunnableOutput]:
        return [MockRunnableOutput("Read quality reports")]

    # Mirrors galaxy-tool-util's output_to_cwl_json, which falls through to a bare
    # `return None` for any collection type it does not know how to translate.
    def mock_output_to_cwl_json(*args, **kwds) -> None:
        return None

    monkeypatch.setattr(activity, "get_outputs", mock_get_outputs)
    monkeypatch.setattr(activity, "output_to_cwl_json", mock_output_to_cwl_json)

    class MockRunnable:
        type = None

    return CollectionRunResponse(
        ctx=cast("PlanemoCliContext", MockContext()),
        runnable=MockRunnable(),
        user_gi=cast("GalaxyInstance", None),
        history_id="a2619fc82a1e31e5",
        log="",
    )


def test_collect_outputs_unsupported_collection_type(unsupported_collection_response, tmp_path):
    """An untranslatable collection type reports the output and type, not a TypeError."""
    with pytest.raises(Exception) as exc_info:
        unsupported_collection_response.collect_outputs(output_directory=str(tmp_path))

    assert not isinstance(exc_info.value, TypeError)
    message = str(exc_info.value)
    assert "Read quality reports" in message
    assert "sample_sheet:paired" in message
