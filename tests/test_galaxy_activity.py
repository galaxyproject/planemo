"""Unit tests for :mod:`planemo.galaxy.activity`."""

import os
from unittest import mock

import pytest
import yaml
from gxformat2.converter import python_to_workflow
from gxformat2.normalize import inputs_normalized
from pydantic import ValidationError

from planemo.galaxy.activity import (
    _execute,
    GalaxyToolRunResponse,
    stage_in,
)
from planemo.runnable import (
    Runnable,
    RunnableType,
)
from .test_utils import (
    create_test_context,
    PROJECT_TEMPLATES_DIR,
    TEST_DATA_DIR,
)

TOOL_RUNNABLE = Runnable(os.path.join(PROJECT_TEMPLATES_DIR, "demo", "cat.xml"), RunnableType.galaxy_tool)
WORKFLOW_RUNNABLE = Runnable(os.path.join(TEST_DATA_DIR, "wf2.ga"), RunnableType.galaxy_workflow)
MULTIPLE_PARAM_WORKFLOW_RUNNABLE = Runnable(
    os.path.join(TEST_DATA_DIR, "wf_multiple_param_inputs.gxwf.yml"), RunnableType.galaxy_workflow
)
COLOR_PARAM_WORKFLOW_RUNNABLE = Runnable(
    os.path.join(TEST_DATA_DIR, "wf_color_param_input.ga"), RunnableType.galaxy_workflow
)


def _gxformat2_normalizes_color_and_directory_uri_inputs():
    steps = {
        str(i): {"id": i, "label": t, "type": "parameter_input", "tool_state": {"parameter_type": t}, "inputs": []}
        for i, t in enumerate(["color", "directory_uri"])
    }
    try:
        inputs_normalized(workflow_dict={"a_galaxy_workflow": "true", "format-version": "0.1", "steps": steps})
    except ValidationError:
        return False
    return True


class ExecutionHalted(Exception):
    """Raised to stop ``_execute`` once the request of interest was captured."""


def _execute_capturing_request(runnable, **kwds):
    """Run ``_execute`` far enough to capture the request it makes of Galaxy."""
    config = mock.MagicMock()
    config.user_gi.tools._post.side_effect = ExecutionHalted()
    config.user_gi.workflows.invoke_workflow.side_effect = ExecutionHalted()
    with mock.patch("planemo.galaxy.activity.stage_in", return_value=({}, "historyid123")):
        with pytest.raises(ExecutionHalted):
            _execute(create_test_context(), config, runnable, job_path=None, **kwds)
    if runnable.type == RunnableType.galaxy_tool:
        return config.user_gi.tools._post.call_args.args[0]
    else:
        return config.user_gi.workflows.invoke_workflow.call_args.kwargs


@pytest.mark.parametrize("use_cache", [True, False])
def test_execute_tool_forwards_use_cache(use_cache):
    payload = _execute_capturing_request(TOOL_RUNNABLE, use_cache=use_cache)
    assert payload["use_cached_job"] is use_cache


@pytest.mark.parametrize("use_cache", [True, False])
def test_execute_workflow_forwards_use_cache(use_cache):
    invoke_kwds = _execute_capturing_request(WORKFLOW_RUNNABLE, use_cache=use_cache)
    assert invoke_kwds["use_cached_job"] is use_cache


def test_execute_does_not_cache_without_use_cache():
    """A caller that leaves use_cache unset gets no caching - e.g. ``planemo test`` by default."""
    assert _execute_capturing_request(TOOL_RUNNABLE)["use_cached_job"] is False
    assert _execute_capturing_request(WORKFLOW_RUNNABLE)["use_cached_job"] is False


def test_execute_tool_no_wait_returns_submission_response():
    config = mock.MagicMock()
    config.log_contents = ""
    api_run_response = {
        "jobs": [{"id": "jobid123"}],
        "outputs": [],
        "output_collections": [],
    }
    config.user_gi.tools._post.return_value = api_run_response

    with mock.patch("planemo.galaxy.activity.stage_in", return_value=({}, "historyid123")):
        with mock.patch("planemo.galaxy.activity._wait_for_job") as wait_for_job:
            response = _execute(create_test_context(), config, TOOL_RUNNABLE, job_path=None, no_wait=True)

    assert isinstance(response, GalaxyToolRunResponse)
    assert response.was_successful
    assert response.job_info is None
    assert response.api_run_response is api_run_response
    assert response.structured_data()["data"]["status"] == "success"
    wait_for_job.assert_not_called()
    config.gi.jobs.show_job.assert_not_called()


def _staging_config():
    """A config whose Galaxy answers dataset uploads and collection creation."""

    def make_post_request(url, payload, files_attached=False):
        if url.endswith("/tools/fetch"):
            return {"outputs": [{"id": "dataset1", "src": "hda"}], "jobs": [{"id": "job1"}]}
        elif url.endswith("/dataset_collections"):
            return {"id": "collection1"}
        raise AssertionError(f"Unexpected POST to {url}")

    config = mock.MagicMock()
    config.version_major = "25.0"
    config.use_path_paste = True
    config.user_gi.url = "http://localhost:8080/api"
    config.user_gi.make_post_request.side_effect = make_post_request
    return config


def _write_multiple_param_job(tmp_path):
    (tmp_path / "a.txt").write_text("a\n")
    job_path = tmp_path / "job.yml"
    job_path.write_text(
        yaml.safe_dump(
            {
                "text_list": ["curated_consensus", "--ex2"],
                "int_list": [1, 2],
                "data_list": [{"class": "File", "path": "a.txt"}],
            }
        )
    )
    return str(job_path)


def _stage_in(runnable, config, job_path):
    with mock.patch("planemo.galaxy.activity._wait_for_job"):
        return stage_in(
            create_test_context(), runnable, config, job_path, history_id="history1", check_uploads_ok=False
        )


EXPECTED_MULTIPLE_PARAM_JOB_DICT = {
    "text_list": ["curated_consensus", "--ex2"],
    "int_list": [1, 2],
    # a list given to a data input is still uploaded as a list collection
    "data_list": {"src": "hdca", "id": "collection1"},
}


def test_stage_in_passes_lists_to_multiple_value_workflow_parameter_inputs(tmp_path):
    job_path = _write_multiple_param_job(tmp_path)
    job_dict, history_id = _stage_in(MULTIPLE_PARAM_WORKFLOW_RUNNABLE, _staging_config(), job_path)
    assert history_id == "history1"
    assert job_dict == EXPECTED_MULTIPLE_PARAM_JOB_DICT


@pytest.mark.parametrize(
    "uri, instance",
    [
        ("gxid://workflows/workflow1", False),
        ("gxid://workflows/workflow1?runnable_path=/does/not/exist.gxwf.yml", False),
        ("gxid://workflow-instance/workflow1", True),
    ],
)
def test_stage_in_passes_lists_to_multiple_value_parameter_inputs_of_remote_workflow(tmp_path, uri, instance):
    """A remote workflow's inputs are read from Galaxy, its runnable_path only locates tests."""
    job_path = _write_multiple_param_job(tmp_path)
    with open(MULTIPLE_PARAM_WORKFLOW_RUNNABLE.path) as f:
        workflow = python_to_workflow(yaml.safe_load(f), workflow_directory=TEST_DATA_DIR)

    def download_workflow(path, params):
        assert (path, params) == ("workflow1/download", {"instance": instance})
        return workflow

    config = _staging_config()
    config.workflow_id_for_runnable.return_value = "workflow1"
    config.user_gi.workflows._get.side_effect = download_workflow
    job_dict, _ = _stage_in(Runnable(uri, RunnableType.galaxy_workflow), config, job_path)
    assert job_dict == EXPECTED_MULTIPLE_PARAM_JOB_DICT


@pytest.mark.parametrize(
    "job, expected",
    [
        ({"color_param": "#ff0000"}, {"color_param": "#ff0000"}),
        ({"data_list": [{"class": "File", "path": "a.txt"}]}, {"data_list": {"src": "hdca", "id": "collection1"}}),
    ],
)
def test_stage_in_without_parameter_lists_skips_workflow_normalization(tmp_path, job, expected):
    """gxformat2 can't normalize color inputs, so only parameter lists may need it."""
    (tmp_path / "a.txt").write_text("a\n")
    job_path = tmp_path / "job.yml"
    job_path.write_text(yaml.safe_dump(job))
    job_dict, _ = _stage_in(COLOR_PARAM_WORKFLOW_RUNNABLE, _staging_config(), str(job_path))
    assert job_dict == expected


@pytest.mark.xfail(
    not _gxformat2_normalizes_color_and_directory_uri_inputs(),
    raises=ValidationError,
    strict=True,
    reason="gxformat2 can't normalize native color and directory_uri inputs (galaxyproject/gxformat2#202)",
)
def test_stage_in_passes_parameter_lists_next_to_color_and_directory_uri_inputs(tmp_path):
    (tmp_path / "a.txt").write_text("a\n")
    job_path = tmp_path / "job.yml"
    params = {"text_list": ["a", "b"], "color_param": "#ff0000", "directory_uri_param": "gxfiles://dropbox/out"}
    job_path.write_text(yaml.safe_dump({**params, "data_list": [{"class": "File", "path": "a.txt"}]}))
    job_dict, _ = _stage_in(COLOR_PARAM_WORKFLOW_RUNNABLE, _staging_config(), str(job_path))
    assert job_dict == {**params, "data_list": {"src": "hdca", "id": "collection1"}}
