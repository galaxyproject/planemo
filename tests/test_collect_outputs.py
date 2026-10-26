"""Tests for output collection in :mod:`planemo.galaxy.activity`."""

import json
from copy import deepcopy
from types import SimpleNamespace
from typing import (
    Any,
    cast,
    Dict,
    List,
    Optional,
    TYPE_CHECKING,
)

import pytest
from galaxy.tool_util.cwl.util import GalaxyOutput

from planemo.galaxy import activity
from planemo.galaxy.activity import GalaxyBaseRunResponse
from planemo.runnable import RunnableType
from planemo.test._check_output import check_output

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


@pytest.mark.parametrize("value", [None, False, 42, "hello", [None, 1], {"answer": 42}])
@pytest.mark.parametrize("nested", [False, True])
@pytest.mark.parametrize("cwl", [False, True])
def test_expression_collection_assertions(monkeypatch, tmp_path, value, nested, cwl):
    """Collect real expression.json content and verify it without losing its file path."""
    content = json.dumps(value)
    dataset = {"id": "expression", "history_content_type": "dataset", "file_ext": "expression.json"}
    leaf = {"element_identifier": "sample", "element_type": "hda", "object": dataset}
    collection = {"id": "inner", "collection_type": "list", "elements": [leaf]}
    element_tests = {"sample": {"asserts": [{"that": "has_text", "text": content}]}}
    if nested:
        collection = {
            "collection_type": "list:list",
            "elements": [{"element_identifier": "group", "element_type": "dataset_collection", "object": collection}],
        }
        element_tests = {"group": {"elements": element_tests}}
    collection.update(id="collection", history_content_type="dataset_collection")

    class ExpressionResponse(CollectionRunResponse):
        def to_galaxy_output(self, runnable_output):
            return GalaxyOutput("history", "dataset_collection", "collection", None)

        def _get_metadata(self, history_content_type, content_id):
            return deepcopy(collection if history_content_type == "dataset_collection" else dataset)

        def download_output_to(self, ctx, details, directory, filename=None):
            path = tmp_path / "expression.json"
            path.write_text(content)
            return str(path)

    runnable = SimpleNamespace(type=RunnableType.cwl_workflow if cwl else RunnableType.galaxy_workflow, path="test.ga")
    response = ExpressionResponse(MockContext(), runnable, None, "history", "")
    monkeypatch.setattr(activity, "get_outputs", lambda *args, **kwds: [MockRunnableOutput("out")])
    response.collect_outputs(output_directory=str(tmp_path))
    output = response.outputs_dict["out"]
    if cwl:
        assert output == ([[value]] if nested else [value])
    else:
        assert check_output(runnable, output, {"name": "out", "element_tests": element_tests}) == []
        failing_tests = {"sample": {"asserts": [{"that": "has_text", "text": "missing-content"}]}}
        if nested:
            failing_tests = {"group": {"elements": failing_tests}}
        problems = check_output(runnable, output, {"name": "out", "element_tests": failing_tests})
        assert len(problems) == 1
        assert "missing-content" in problems[0]


@pytest.mark.parametrize("value", [None, False, 42, "hello", [None, 1], {"answer": 42}])
def test_collection_element_without_file_properties(value):
    """A malformed element reports an output problem instead of aborting the test run."""
    output = {
        "collection_type": "list",
        "elements": [{"element_identifier": "sample", "element_type": "hda", "object": {}, "_output_object": value}],
    }
    tests = {"name": "out", "element_tests": {"sample": {"asserts": [{"that": "has_text", "text": "null"}]}}}
    problems = check_output(SimpleNamespace(path="test.ga"), output, tests)
    assert problems == ["Expected file properties for element 'sample' of output collection 'out'"]


def test_mixed_collection_assertions(monkeypatch, tmp_path):
    """Different expression files and an ordinary dataset retain their own content."""
    contents = {"skipped": "null", "number": "42", "ordinary": "plain text"}
    datasets = {
        name: {
            "id": name,
            "history_content_type": "dataset",
            "file_ext": "txt" if name == "ordinary" else "expression.json",
        }
        for name in contents
    }
    metadata = {
        "id": "collection",
        "history_content_type": "dataset_collection",
        "collection_type": "list",
        "elements": [
            {"element_identifier": name, "element_type": "hda", "object": dataset} for name, dataset in datasets.items()
        ],
    }

    class MixedResponse(CollectionRunResponse):
        def to_galaxy_output(self, runnable_output):
            return GalaxyOutput("history", "dataset_collection", "collection", None)

        def _get_metadata(self, history_content_type, content_id):
            return deepcopy(metadata if history_content_type == "dataset_collection" else datasets[content_id])

        def _get_extra_files(self, dataset_details):
            return []

        def download_output_to(self, ctx, details, directory, filename=None):
            path = tmp_path / details["id"]
            path.write_text(contents[details["id"]])
            return str(path)

    runnable = SimpleNamespace(type=RunnableType.galaxy_workflow, path="test.ga")
    response = MixedResponse(MockContext(), runnable, None, "history", "")
    monkeypatch.setattr(activity, "get_outputs", lambda *args, **kwds: [MockRunnableOutput("out")])
    output = response.collect_outputs(output_directory=str(tmp_path), output_id="out")
    tests = {
        "name": "out",
        "element_tests": {
            name: {"asserts": [{"that": "has_text", "text": content}]} for name, content in contents.items()
        },
    }
    assert check_output(runnable, output, tests) == []
    tests["element_tests"]["number"]["asserts"][0]["text"] = "null"
    assert len(check_output(runnable, output, tests)) == 1
