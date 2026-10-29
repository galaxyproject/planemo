"""Unit tests for engines and runnables."""

import os

import bioblend

from planemo.engine import engine_context
from planemo.engine.galaxy import log_service_logs_on_failure
from planemo.engine.interface import BaseEngine
from planemo.runnable import (
    for_path,
    get_outputs,
    RunnableType,
    SuccessfulRunResponse,
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


class _ShutDownRunResponse(SuccessfulRunResponse):
    """Fetches every output from its engine and can't once that has shut down - like a
    Galaxy run response whose output cache misses after Galaxy has stopped."""

    def __init__(self, runnable, outputs):
        super().__init__(runnable)
        self._outputs = outputs
        self.engine_running = True

    @property
    def log(self):
        return None

    @property
    def job_info(self):
        return None

    @property
    def invocation_details(self):
        return None

    @property
    def outputs_dict(self):
        return {}

    def get_output(self, output_id):
        if not self.engine_running:
            raise bioblend.ConnectionError("Connection refused")
        return self._outputs.get(output_id)


class _ShutDownEngine(BaseEngine):
    """Calls the output collectors while running, like GalaxyEngine, then shuts down."""

    handled_runnable_types = [RunnableType.galaxy_workflow]

    def __init__(self, ctx, outputs, **kwds):
        super().__init__(ctx, **kwds)
        self._outputs = outputs

    def _run(self, runnables, job_paths, output_collectors=None, test_timeout=None):
        run_responses = [_ShutDownRunResponse(runnable, self._outputs) for runnable in runnables]
        for run_response, collect_output in zip(run_responses, output_collectors or []):
            collect_output(run_response)
        for run_response in run_responses:
            run_response.engine_running = False
        return run_responses


def test_test_results_do_not_need_engine_after_shutdown(tmp_path):
    """Structured test data is built while the engine runs, not rebuilt after it has shut down."""
    output_path = tmp_path / "wf_output_1.txt"
    output_path.write_text("Hello World!\n")
    runnable = for_path(os.path.join(TEST_DATA_DIR, "wf20_undefined_output.yml"))
    engine = _ShutDownEngine(create_test_context(), outputs={"wf_output_1": {"path": str(output_path)}})

    structured_results = engine.test([runnable], test_timeout=None)

    (test_data,) = structured_results.structured_data["tests"]
    assert test_data["data"]["status"] == "failure"
    assert test_data["data"]["output_problems"] == ["Expected output [wf_output_typo] not found in results."]
