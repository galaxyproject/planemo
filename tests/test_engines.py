"""Unit tests for engines and runnables."""

import contextlib
import os
from types import SimpleNamespace
from unittest.mock import (
    Mock,
    patch,
)

from planemo.engine import engine_context
from planemo.engine.galaxy import (
    DockerizedManagedGalaxyEngine,
    InstalledGalaxyEngine,
    LocalManagedGalaxyEngine,
    log_service_logs_on_failure,
)
from planemo.engine.interface import BaseEngine
from planemo.runnable import (
    for_path,
    get_outputs,
    RunnableType,
)
from planemo.test.results import StructuredData
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


class _RecordingManagedGalaxyEngine(LocalManagedGalaxyEngine):
    def __init__(self, ctx, config):
        super().__init__(ctx)
        self.config = config
        self.served = []

    @contextlib.contextmanager
    def ensure_runnables_served(self, runnables):
        self.served.append(list(runnables))
        yield self.config


def test_managed_galaxy_engines_can_serve_test_results():
    assert LocalManagedGalaxyEngine.can_serve_test_results
    assert InstalledGalaxyEngine.can_serve_test_results
    assert DockerizedManagedGalaxyEngine.can_serve_test_results


def test_managed_test_context_reuses_one_server_until_reporting_finishes():
    config = SimpleNamespace(galaxy_url="http://127.0.0.1:9090")
    engine = _RecordingManagedGalaxyEngine(create_test_context(), config)
    result = StructuredData(data={"version": "0.1", "tests": []})
    result.calculate_summary_data()

    with patch.object(BaseEngine, "test", return_value=result):
        with engine.test_context(["one", "two"], test_timeout=17, keep_alive=True) as actual_result:
            assert actual_result is result
            with engine._served_config(["one"]) as active_config:
                assert active_config is config

    assert engine.served == [["one", "two"]]


def test_managed_test_server_stops_cleanly_after_interrupt():
    config = SimpleNamespace(galaxy_url="http://127.0.0.1:9090")
    engine = _RecordingManagedGalaxyEngine(create_test_context(), config)

    with (
        patch.object(BaseEngine, "test", return_value=Mock()),
        patch("planemo.engine.galaxy.sleep_for_serve", side_effect=KeyboardInterrupt) as sleep,
        engine.test_context([], test_timeout=17, keep_alive=True),
    ):
        engine.serve_test_results()

    sleep.assert_called_once_with()


def test_embedded_tool_test_histories_are_only_preserved_while_serving():
    config = SimpleNamespace(
        galaxy_url="http://127.0.0.1:9090",
        master_api_key="master-key",
        user_api_key="user-key",
        service_log_contents={},
    )
    engine = _RecordingManagedGalaxyEngine(create_test_context(), config)
    test_case = SimpleNamespace(tool_id="cat", test_index=0, tool_version="1.0")

    with patch("planemo.engine.galaxy.interactor.verify_tool") as verify_tool:
        engine._run_galaxy_tool_test_case(config, test_case, 17, Mock())
        with (
            patch.object(BaseEngine, "test", return_value=Mock()),
            engine.test_context([], test_timeout=17, keep_alive=True),
        ):
            engine._run_galaxy_tool_test_case(config, test_case, 17, Mock())

    assert verify_tool.call_args_list[0].kwargs["no_history_cleanup"] is False
    assert verify_tool.call_args_list[1].kwargs["no_history_cleanup"] is True
