"""Unit tests for engines and runnables."""

import os

import pytest
import responses

from planemo.engine import engine_context
from planemo.engine.galaxy import (
    log_service_logs_on_failure,
    PlanemoGalaxyInteractorApi,
)
from planemo.runnable import (
    for_path,
    get_outputs,
    RunnableType,
)
from .test_utils import (
    create_test_context,
    TEST_DATA_DIR,
)

A_CWL_TOOL = os.path.join(TEST_DATA_DIR, "tools", "ok-cat1-tool.cwl")
A_CWL_WORKFLOW = os.path.join(TEST_DATA_DIR, "count-lines2-wf.cwl")

A_GALAXY_TOOL = os.path.join(TEST_DATA_DIR, "tools", "ok_select_param.xml")
A_GALAXY_GA_WORKFLOW = os.path.join(TEST_DATA_DIR, "test_workflow_1.ga")
A_GALAXY_YAML_WORKFLOW = os.path.join(TEST_DATA_DIR, "wf1.gxwf.yml")

CAN_HANDLE = {
    "galaxy": {
        A_CWL_TOOL: True,
        A_CWL_WORKFLOW: True,
        A_GALAXY_TOOL: True,
        A_GALAXY_GA_WORKFLOW: True,
        A_GALAXY_YAML_WORKFLOW: True,
    },
    "cwltool": {
        A_CWL_TOOL: True,
        A_CWL_WORKFLOW: True,
        A_GALAXY_TOOL: False,
        A_GALAXY_GA_WORKFLOW: False,
        A_GALAXY_YAML_WORKFLOW: False,
    },
}


def test_can_handle():
    ctx = create_test_context()
    for engine_type in ["galaxy", "cwltool"]:
        with engine_context(ctx, engine=engine_type) as e:
            for key, value in CAN_HANDLE[engine_type].items():
                assert bool(e.can_run(for_path(key))) is value


def test_outputs():
    outputs = get_outputs(for_path(A_CWL_WORKFLOW))
    assert len(outputs) == 1
    output_id = outputs[0].get_id()
    assert output_id == "count_output"


def test_runnable_types():
    assert RunnableType.galaxy_tool.is_galaxy_artifact
    assert RunnableType.galaxy_workflow.is_galaxy_artifact
    assert not RunnableType.cwl_tool.is_galaxy_artifact
    assert not RunnableType.cwl_workflow.is_galaxy_artifact


class _RecordingContext:
    def __init__(self):
        self.messages = []

    def log(self, msg, *args):
        self.messages.append(msg)


class _ConfigWithServiceLogs:
    service_log_contents = {"celery.log": "Task galaxy.fetch_data raised unexpected"}


def test_service_logs_not_logged_when_tests_pass():
    ctx = _RecordingContext()
    log_service_logs_on_failure(ctx, _ConfigWithServiceLogs(), [{"data": {"status": "success"}}])
    assert ctx.messages == []


def test_service_logs_logged_when_test_fails():
    ctx = _RecordingContext()
    log_service_logs_on_failure(ctx, _ConfigWithServiceLogs(), [{"data": {"status": "error"}}])
    assert len(ctx.messages) == 1
    assert "celery.log" in ctx.messages[0]
    assert "raised unexpected" in ctx.messages[0]


def test_service_logs_logged_when_no_result_registered():
    """verify_tool blew up before registering anything - still want the logs."""
    ctx = _RecordingContext()
    log_service_logs_on_failure(ctx, _ConfigWithServiceLogs(), [])
    assert len(ctx.messages) == 1


@responses.activate
def test_embedded_test_data_download_with_user_key_only():
    galaxy_interactor = PlanemoGalaxyInteractorApi(
        galaxy_url="https://example.org", master_api_key=None, api_key="user-key"
    )
    responses.get("https://example.org/api/version", json={"version_major": "25.1"})
    responses.get("https://example.org/api/tools/cat/test_data_download?filename=input.txt", body="test input")

    assert galaxy_interactor.test_data_download("cat", "input.txt") == b"test input"
    assert all(call.request.headers["x-api-key"] == "user-key" for call in responses.calls)


@responses.activate
def test_embedded_test_server_path_requires_admin_key():
    galaxy_interactor = PlanemoGalaxyInteractorApi(
        galaxy_url="https://example.org", master_api_key=None, api_key="user-key"
    )
    responses.get("https://example.org/api/version", json={"version_major": "18.09"})

    with pytest.raises(RuntimeError, match="--galaxy_admin_key") as exc:
        galaxy_interactor.test_data_download("cat", "input.txt")

    assert "server-side test data paths" in str(exc.value)
    assert "user-key" not in str(exc.value)
    assert len(responses.calls) == 1


@responses.activate
def test_embedded_test_server_path_preserves_admin_credentials():
    galaxy_interactor = PlanemoGalaxyInteractorApi(
        galaxy_url="https://example.org", master_api_key="admin-key", api_key="user-key"
    )
    responses.get(
        "https://example.org/api/tools/cat/test_data_path?filename=input.txt&tool_version=1.0",
        json="/galaxy/test-data/input.txt",
    )

    assert galaxy_interactor.test_data_path("cat", "input.txt", tool_version="1.0") == "/galaxy/test-data/input.txt"
    assert responses.calls[0].request.headers["x-api-key"] == "admin-key"
