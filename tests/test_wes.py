"""Tests for the GA4GH WES client and ``planemo run --wes`` plumbing."""

import json
from pathlib import Path
from unittest import mock

import pytest
import responses
from click.testing import CliRunner

from planemo.cli import PlanemoCliContext
from planemo.commands.cmd_run import cli as run_cli
from planemo.galaxy.activity import (
    _is_unstaged_file_input,
    _wait_for_wes_run,
    _wes_execute,
    _wes_workflow_params,
    _wes_workflow_type,
    GalaxyWorkflowRunResponse,
    WesRunResponse,
)
from planemo.galaxy.wes import (
    detect_workflow_type,
    FAILURE_STATES,
    is_failure,
    is_success,
    is_terminal,
    STATE_COMPLETE,
    TERMINAL_STATES,
    WesClient,
    WesError,
)
from planemo.runnable import (
    Runnable,
    RunnableType,
)

GALAXY_URL = "http://localhost:8080"
WES_RUNS = f"{GALAXY_URL}/ga4gh/wes/v1/runs"
DATA_DIR = Path(__file__).parent / "data"
GA_WORKFLOW = str(DATA_DIR / "test_workflow_1.ga")
CWL_WORKFLOW = str(DATA_DIR / "count-lines2-wf.cwl")
TOOL_XML = str(DATA_DIR / "cat_list.xml")


class TestDetectWorkflowType:
    def test_format2_yaml(self):
        assert detect_workflow_type("class: GalaxyWorkflow\nsteps: {}") == "gx_workflow_format2"

    def test_format2_json(self):
        assert detect_workflow_type(json.dumps({"class": "GalaxyWorkflow"})) == "gx_workflow_format2"

    def test_native_ga_json(self):
        assert detect_workflow_type(json.dumps({"a_galaxy_workflow": "true", "steps": {}})) == "gx_workflow_ga"

    def test_native_ga_default(self):
        assert detect_workflow_type("steps:\n  0: {}") == "gx_workflow_ga"


class TestWesClient:
    @responses.activate
    def test_submit_run_sends_params_and_returns_run_id(self):
        responses.add(responses.POST, WES_RUNS, json={"run_id": "abc123"}, status=200)
        client = WesClient(GALAXY_URL, api_key="key123")
        result = client.submit_run(
            workflow_type="gx_workflow_ga",
            workflow_url="gxworkflow://abc",
            params={"input1": 42},
            engine_parameters={"history_name": "WES Run"},
        )
        assert result["run_id"] == "abc123"
        request = responses.calls[0].request
        assert request.headers["x-api-key"] == "key123"
        body = request.body
        assert "workflow_type=gx_workflow_ga" in body
        # JSON-encoded params are url-encoded into the form body.
        assert "input1" in body

    @responses.activate
    def test_submit_run_requires_workflow_url(self):
        client = WesClient(GALAXY_URL)
        with pytest.raises(ValueError):
            client.submit_run(workflow_type="gx_workflow_ga", workflow_url="")

    @responses.activate
    def test_get_run_status(self):
        responses.add(
            responses.GET,
            f"{WES_RUNS}/abc123/status",
            json={"run_id": "abc123", "state": "RUNNING"},
            status=200,
        )
        client = WesClient(GALAXY_URL, api_key="key123")
        assert client.get_run_status("abc123")["state"] == "RUNNING"

    @responses.activate
    def test_error_raises_wes_error_with_status_code(self):
        responses.add(responses.GET, f"{WES_RUNS}/missing/status", json={"err_msg": "nope"}, status=404)
        client = WesClient(GALAXY_URL, api_key="key123")
        with pytest.raises(WesError) as exc_info:
            client.get_run_status("missing")
        assert exc_info.value.status_code == 404
        assert "nope" in str(exc_info.value)


class TestStateConstants:
    def test_complete_is_terminal_not_failure(self):
        assert STATE_COMPLETE in TERMINAL_STATES
        assert STATE_COMPLETE not in FAILURE_STATES

    def test_failures_are_terminal(self):
        assert FAILURE_STATES <= TERMINAL_STATES


class TestWorkflowParams:
    def test_scalar_params_allowed(self, tmp_path):
        job = tmp_path / "job.yml"
        job.write_text("input1: 5\ninput2: hello\n")
        assert _wes_workflow_params(str(job)) == {"input1": 5, "input2": "hello"}

    def test_prestaged_hda_reference_allowed(self, tmp_path):
        job = tmp_path / "job.yml"
        job.write_text("input1:\n  src: hda\n  id: abc123\n")
        params = _wes_workflow_params(str(job))
        assert params["input1"] == {"src": "hda", "id": "abc123"}

    def test_local_file_input_rejected(self, tmp_path):
        job = tmp_path / "job.yml"
        job.write_text("input1:\n  class: File\n  path: ./data.txt\n")
        with pytest.raises(Exception, match="data-staging"):
            _wes_workflow_params(str(job))

    def test_location_input_rejected(self, tmp_path):
        job = tmp_path / "job.yml"
        job.write_text("input1:\n  location: http://example.com/data.txt\n")
        with pytest.raises(Exception, match="data-staging"):
            _wes_workflow_params(str(job))

    def test_is_unstaged_detects_nested_list(self):
        assert _is_unstaged_file_input([{"class": "File", "path": "x"}])
        assert not _is_unstaged_file_input([{"src": "hda", "id": "1"}])
        assert not _is_unstaged_file_input("scalar")


class TestStateHelpers:
    def test_is_success(self):
        assert is_success("COMPLETE")
        assert not is_success("RUNNING")
        assert not is_success("SYSTEM_ERROR")

    def test_is_failure(self):
        assert is_failure("SYSTEM_ERROR")
        assert is_failure("EXECUTOR_ERROR")
        assert not is_failure("COMPLETE")
        assert not is_failure("RUNNING")

    def test_is_terminal(self):
        assert is_terminal("COMPLETE")
        assert is_terminal("CANCELED")
        assert not is_terminal("RUNNING")


class TestWesWorkflowType:
    def test_native_ga_path(self):
        runnable = Runnable(GA_WORKFLOW, RunnableType.galaxy_workflow)
        assert _wes_workflow_type(runnable) == "gx_workflow_ga"

    def test_format2_path(self, tmp_path):
        wf = tmp_path / "wf.gxwf.yml"
        wf.write_text("class: GalaxyWorkflow\nsteps: {}\n")
        runnable = Runnable(str(wf), RunnableType.galaxy_workflow)
        assert _wes_workflow_type(runnable) == "gx_workflow_format2"

    def test_remote_instance_uri_without_path_defaults_to_ga(self):
        runnable = Runnable("gxid://workflow-instance/abc123", RunnableType.galaxy_workflow)
        assert not runnable.has_path
        assert _wes_workflow_type(runnable) == "gx_workflow_ga"


class TestWaitForWesRun:
    def test_returns_terminal_state(self):
        client = mock.Mock()
        client.get_run_status.side_effect = [{"state": "RUNNING"}, {"state": "COMPLETE"}]
        with mock.patch("planemo.galaxy.activity.time.sleep"):
            state = _wait_for_wes_run(mock.MagicMock(), client, "run1")
        assert state == "COMPLETE"
        assert client.get_run_status.call_count == 2

    def test_times_out(self):
        client = mock.Mock()
        client.get_run_status.return_value = {"state": "RUNNING"}
        with mock.patch("planemo.galaxy.activity.time.sleep"):
            with pytest.raises(Exception, match="Timed out"):
                _wait_for_wes_run(mock.MagicMock(), client, "run1", test_timeout=0)


class TestWesRunResponseSuccess:
    def _response(self, wes_state, no_wait=False):
        with mock.patch.object(GalaxyWorkflowRunResponse, "collect_invocation_details", return_value={}):
            return WesRunResponse(
                ctx=mock.MagicMock(),
                runnable=Runnable(GA_WORKFLOW, RunnableType.galaxy_workflow),
                user_gi=mock.MagicMock(),
                history_id="hist1",
                workflow_id="wf1",
                invocation_id="inv1",
                wes_state=wes_state,
                log="log",
                no_wait=no_wait,
            )

    def test_complete_is_successful(self):
        assert self._response("COMPLETE").was_successful

    def test_failure_is_not_successful(self):
        assert not self._response("SYSTEM_ERROR").was_successful

    def test_running_is_not_successful_when_waited(self):
        assert not self._response("RUNNING").was_successful

    def test_no_wait_is_successful_regardless(self):
        assert self._response("RUNNING", no_wait=True).was_successful


def _fake_config(workflow_ref_id="wfid", history_id="hist1"):
    config = mock.MagicMock()
    config.galaxy_url = GALAXY_URL
    config.user_api_key = "key123"
    config.workflow_id_for_runnable.return_value = workflow_ref_id
    config.user_gi.invocations.show_invocation.return_value = {
        "history_id": history_id,
        "workflow_id": "wfinstance",
    }
    return config


class TestWesExecute:
    @mock.patch("planemo.galaxy.activity.summarize_history")
    @mock.patch("planemo.galaxy.activity.WesRunResponse")
    @mock.patch("planemo.galaxy.activity.WesClient")
    def test_success_submits_via_gxworkflow_uri(self, MockClient, MockResponse, mock_summarize, tmp_path):
        job = tmp_path / "job.yml"
        job.write_text("param1: 7\n")
        client = MockClient.return_value
        client.submit_run.return_value = {"run_id": "inv1"}
        client.get_run_status.return_value = {"state": "COMPLETE"}
        runnable = Runnable(GA_WORKFLOW, RunnableType.galaxy_workflow)

        result = _wes_execute(mock.MagicMock(), _fake_config(), runnable, str(job))

        submit_kwargs = client.submit_run.call_args.kwargs
        assert submit_kwargs["workflow_url"] == "gxworkflow://wfid"
        assert submit_kwargs["workflow_type"] == "gx_workflow_ga"
        assert submit_kwargs["params"] == {"param1": 7}
        response_kwargs = MockResponse.call_args.kwargs
        assert response_kwargs["wes_state"] == "COMPLETE"
        assert response_kwargs["invocation_id"] == "inv1"
        assert response_kwargs["history_id"] == "hist1"
        mock_summarize.assert_not_called()
        assert result is MockResponse.return_value

    @mock.patch("planemo.galaxy.activity.summarize_history")
    @mock.patch("planemo.galaxy.activity.WesRunResponse")
    @mock.patch("planemo.galaxy.activity.WesClient")
    def test_failure_summarizes_history(self, MockClient, MockResponse, mock_summarize, tmp_path):
        job = tmp_path / "job.yml"
        job.write_text("param1: 7\n")
        client = MockClient.return_value
        client.submit_run.return_value = {"run_id": "inv1"}
        client.get_run_status.return_value = {"state": "SYSTEM_ERROR"}
        runnable = Runnable(GA_WORKFLOW, RunnableType.galaxy_workflow)

        _wes_execute(mock.MagicMock(), _fake_config(), runnable, str(job))

        mock_summarize.assert_called_once()
        assert MockResponse.call_args.kwargs["wes_state"] == "SYSTEM_ERROR"

    @mock.patch("planemo.galaxy.activity.summarize_history")
    @mock.patch("planemo.galaxy.activity.WesRunResponse")
    @mock.patch("planemo.galaxy.activity.WesClient")
    def test_no_wait_does_not_poll_or_summarize(self, MockClient, MockResponse, mock_summarize, tmp_path):
        job = tmp_path / "job.yml"
        job.write_text("param1: 7\n")
        client = MockClient.return_value
        client.submit_run.return_value = {"run_id": "inv1"}
        client.get_run_status.return_value = {"state": "QUEUED"}
        runnable = Runnable(GA_WORKFLOW, RunnableType.galaxy_workflow)

        _wes_execute(mock.MagicMock(), _fake_config(), runnable, str(job), no_wait=True)

        client.get_run_status.assert_called_once()
        mock_summarize.assert_not_called()
        assert MockResponse.call_args.kwargs["no_wait"] is True

    @mock.patch("planemo.galaxy.activity.summarize_history")
    @mock.patch("planemo.galaxy.activity.WesRunResponse")
    @mock.patch("planemo.galaxy.activity.WesClient")
    def test_instance_uri_appends_instance_flag(self, MockClient, MockResponse, mock_summarize, tmp_path):
        job = tmp_path / "job.yml"
        job.write_text("param1: 7\n")
        client = MockClient.return_value
        client.submit_run.return_value = {"run_id": "inv1"}
        client.get_run_status.return_value = {"state": "COMPLETE"}
        runnable = Runnable("gxid://workflow-instance/abc", RunnableType.galaxy_workflow)

        _wes_execute(mock.MagicMock(), _fake_config(workflow_ref_id="abc"), runnable, str(job))

        assert client.submit_run.call_args.kwargs["workflow_url"] == "gxworkflow://abc?instance=true"

    @mock.patch("planemo.galaxy.activity.summarize_history")
    @mock.patch("planemo.galaxy.activity.WesClient")
    def test_download_outputs_uses_galaxy_api(self, MockClient, mock_summarize, tmp_path):
        job = tmp_path / "job.yml"
        job.write_text("param1: 7\n")
        out_dir = tmp_path / "outs"
        client = MockClient.return_value
        client.submit_run.return_value = {"run_id": "inv1"}
        client.get_run_status.return_value = {"state": "COMPLETE"}
        runnable = Runnable(GA_WORKFLOW, RunnableType.galaxy_workflow)

        with mock.patch("planemo.galaxy.activity.WesRunResponse") as MockResponse:
            _wes_execute(
                mock.MagicMock(),
                _fake_config(),
                runnable,
                str(job),
                download_outputs=True,
                output_directory=str(out_dir),
            )
        MockResponse.return_value.collect_outputs.assert_called_once_with(str(out_dir))


class TestRunCommandOption:
    def test_wes_option_registered(self):
        option_names = {opt.name for opt in run_cli.params}
        assert "wes" in option_names

    def _invoke(self, args):
        return CliRunner().invoke(run_cli, args, obj=PlanemoCliContext(), standalone_mode=False)

    def test_wes_rejects_cwl_workflow(self):
        result = self._invoke([CWL_WORKFLOW, CWL_WORKFLOW, "--wes"])
        assert result.exit_code != 0
        assert "Galaxy workflows" in str(result.exception)

    def test_wes_rejects_tool(self):
        result = self._invoke([TOOL_XML, TOOL_XML, "--wes"])
        assert result.exit_code != 0
        assert "Galaxy workflows" in str(result.exception)

    def test_wes_rejects_non_galaxy_engine(self):
        result = self._invoke([GA_WORKFLOW, GA_WORKFLOW, "--wes", "--engine", "toil"])
        assert result.exit_code != 0
        assert "external_galaxy" in str(result.exception)
