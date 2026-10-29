"""Unit tests for :mod:`planemo.galaxy.activity`."""

import os
from unittest import mock

import pytest

from planemo.galaxy.activity import (
    _execute,
    GalaxyToolRunResponse,
    INFORMAL_PARAMETER_SENTINEL,
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


class ExecutionHalted(Exception):
    """Raised to stop ``_execute`` once the request of interest was captured."""


RUN_DATA = {
    "steps": [
        {"step_type": "parameter_input", "step_label": "declared", "replacement_parameters": []},
        {"step_type": "tool", "step_label": "cat", "replacement_parameters": ["undeclared", "declared"]},
        {"step_type": "subworkflow", "step_label": "sub", "replacement_parameters": ["nested"]},
    ]
}


def _execute_capturing_request(runnable, run_data=None, **kwds):
    """Run ``_execute`` far enough to capture the request it makes of Galaxy."""
    config = mock.MagicMock()
    config.user_gi.workflows._get.return_value = run_data or {"steps": []}
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


def test_execute_workflow_fills_undeclared_parameters_with_sentinel():
    with mock.patch("planemo.galaxy.activity.warn") as warn:
        invoke_kwds = _execute_capturing_request(WORKFLOW_RUNNABLE, run_data=RUN_DATA)
    assert invoke_kwds["replacement_params"] == {
        "nested": INFORMAL_PARAMETER_SENTINEL,
        "undeclared": INFORMAL_PARAMETER_SENTINEL,
    }
    warned = [call.args[0] for call in warn.call_args_list]
    assert warned == [
        f'Workflow uses undeclared parameter ${{{name}}}; sending "{INFORMAL_PARAMETER_SENTINEL}". '
        f"Declare a parameter input labelled {name}."
        for name in ("nested", "undeclared")
    ]


def test_execute_workflow_without_undeclared_parameters_sends_none():
    assert _execute_capturing_request(WORKFLOW_RUNNABLE)["replacement_params"] is None


def test_execute_workflow_invokes_without_params_when_run_data_fails():
    config = mock.MagicMock()
    config.user_gi.workflows._get.side_effect = Exception("Following tools missing")
    config.user_gi.workflows.invoke_workflow.side_effect = ExecutionHalted()
    with mock.patch("planemo.galaxy.activity.stage_in", return_value=({}, "historyid123")):
        with mock.patch("planemo.galaxy.activity.warn") as warn:
            with pytest.raises(ExecutionHalted):
                _execute(create_test_context(), config, WORKFLOW_RUNNABLE, job_path=None)
    assert config.user_gi.workflows.invoke_workflow.call_args.kwargs["replacement_params"] is None
    assert "Following tools missing" in warn.call_args.args[0]


def test_execute_workflow_rejects_job_supplying_undeclared_parameter():
    config = mock.MagicMock()
    config.user_gi.workflows._get.return_value = RUN_DATA
    with mock.patch("planemo.galaxy.activity.stage_in", return_value=({"undeclared": "x"}, "historyid123")):
        with pytest.raises(Exception, match=r"undeclared workflow parameter\(s\) undeclared;"):
            _execute(create_test_context(), config, WORKFLOW_RUNNABLE, job_path=None)
    config.user_gi.workflows.invoke_workflow.assert_not_called()
