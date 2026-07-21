"""Unit tests for :mod:`planemo.galaxy.activity`."""

import os
from unittest import mock

import pytest

from planemo.galaxy.activity import _execute
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
    """Callers other than ``planemo run`` (e.g. ``planemo test``) never set use_cache."""
    assert _execute_capturing_request(TOOL_RUNNABLE)["use_cached_job"] is False
    assert _execute_capturing_request(WORKFLOW_RUNNABLE)["use_cached_job"] is False
