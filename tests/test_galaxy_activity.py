"""Unit tests for :mod:`planemo.galaxy.activity`."""

import os
from unittest import mock

import pytest
import responses

from planemo.galaxy.activity import (
    _execute,
    GalaxyToolRunResponse,
)
from planemo.galaxy.config import external_galaxy_config
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


@pytest.mark.parametrize("admin_key", [None, "admin-key"])
@responses.activate
def test_execute_tool_fetches_job_details_with_available_credentials(admin_key):
    ctx = create_test_context()
    responses.post("https://example.org/api/tools", json={"jobs": [{"id": "jobid123"}], "outputs": []})
    responses.get("https://example.org/api/jobs/jobid123", json={"state": "ok", "stdout": "tool output"})

    with external_galaxy_config(
        ctx, [], galaxy_url="https://example.org", galaxy_user_key="user-key", galaxy_admin_key=admin_key
    ) as config:
        with (
            mock.patch("planemo.galaxy.activity.stage_in", return_value=({}, "historyid123")),
            mock.patch("planemo.galaxy.activity._verified_tool_id", return_value="cat"),
            mock.patch("planemo.galaxy.activity._wait_for_job", return_value="ok"),
        ):
            result = _execute(ctx, config, TOOL_RUNNABLE, job_path=None)

    assert result.was_successful
    assert result.job_info["stdout"] == "tool output"
    assert responses.calls[0].request.headers["x-api-key"] == "user-key"
    assert responses.calls[1].request.headers.get("x-api-key") == (admin_key or "user-key")
