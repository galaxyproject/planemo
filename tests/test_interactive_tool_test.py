"""Unit tests for interactive tool testing (``planemo.galaxy.test.interactive``)."""

import os
from unittest import mock

import pytest

from planemo.galaxy.test.interactive import (
    _split_proxied_url,
    inputs_from_test_case,
    is_interactive_tool,
    resolve_interactive_inputs,
    run_interactive_tool_test,
)
from planemo.galaxy.test.interactive import test_case_count as count_test_cases
from planemo.galaxy.test.interactive import tool_id_for_path
from planemo.test.models import PlanemoTestReport

TESTS_DIR = os.path.dirname(__file__)
INTERACTIVE_TOOL = os.path.join(TESTS_DIR, "data", "tools", "interactivetool_simple.xml")
INTERACTIVE_TOOL_WITH_INPUTS = os.path.join(TESTS_DIR, "data", "tools", "interactivetool_with_inputs.xml")
REGULAR_TOOL = os.path.join(TESTS_DIR, "data", "tools", "bwa_without_requirements.xml")


def test_is_interactive_tool_detects_interactive():
    assert is_interactive_tool(INTERACTIVE_TOOL) is True


def test_is_interactive_tool_detects_interactive_with_inputs():
    assert is_interactive_tool(INTERACTIVE_TOOL_WITH_INPUTS) is True


def test_is_interactive_tool_false_for_regular():
    assert is_interactive_tool(REGULAR_TOOL) is False


def test_is_interactive_tool_false_for_missing_file():
    assert is_interactive_tool(os.path.join(TESTS_DIR, "data", "tools", "does_not_exist.xml")) is False


def test_tool_id_for_path():
    assert tool_id_for_path(INTERACTIVE_TOOL) == "interactivetool_simple"


def test_test_case_count():
    assert count_test_cases(INTERACTIVE_TOOL_WITH_INPUTS) == 2
    # The simple tool declares no <tests>.
    assert count_test_cases(INTERACTIVE_TOOL) == 0


def test_inputs_from_test_case_valid_indexes():
    assert inputs_from_test_case(INTERACTIVE_TOOL_WITH_INPUTS, 1) == {"greeting": "moo cow"}
    assert inputs_from_test_case(INTERACTIVE_TOOL_WITH_INPUTS, 2) == {"greeting": "baa sheep"}


def test_inputs_from_test_case_out_of_range():
    with pytest.raises(ValueError, match="out of range"):
        inputs_from_test_case(INTERACTIVE_TOOL_WITH_INPUTS, 3)


def test_inputs_from_test_case_no_tests():
    with pytest.raises(ValueError, match="0 test case"):
        inputs_from_test_case(INTERACTIVE_TOOL, 1)


def test_resolve_interactive_inputs_test_index_only():
    assert resolve_interactive_inputs(INTERACTIVE_TOOL_WITH_INPUTS, 1, None) == {"greeting": "moo cow"}


def test_resolve_interactive_inputs_manual_overrides_test_index():
    assert resolve_interactive_inputs(INTERACTIVE_TOOL_WITH_INPUTS, 1, {"greeting": "override"}) == {
        "greeting": "override"
    }


def test_resolve_interactive_inputs_manual_only():
    assert resolve_interactive_inputs(INTERACTIVE_TOOL_WITH_INPUTS, None, {"greeting": "x"}) == {"greeting": "x"}


def test_resolve_interactive_inputs_none():
    assert resolve_interactive_inputs(INTERACTIVE_TOOL_WITH_INPUTS, None, None) == {}


@pytest.mark.parametrize(
    "target_url, expected",
    (
        (
            "http://k-t.ep.interactivetool.localhost:8090/",
            ("http://localhost:8090", "k-t.ep.interactivetool.localhost:8090"),
        ),
        (
            "http://k-t.ep.interactivetool.0.0.0.0:8090/foo/bar",
            ("http://0.0.0.0:8090", "k-t.ep.interactivetool.0.0.0.0:8090"),
        ),
        ("http://localhost:8090/plain", ("http://localhost:8090/plain", None)),
    ),
)
def test_split_proxied_url(target_url, expected):
    assert _split_proxied_url(target_url) == expected


class _FakeResponse:
    def __init__(self, payload, status=200):
        self._payload = payload
        self.status_code = status

    def json(self):
        return self._payload

    def raise_for_status(self):
        if self.status_code >= 400:
            raise Exception(f"HTTP {self.status_code}")


def _make_config(entry_points=None, job_state_value="running", target=None, launch_error=None):
    """Build a mock config whose ``user_gi`` behaves like the GxIT API."""
    if entry_points is None:
        entry_points = [{"id": "ep1", "name": "Simple IT", "active": True}]
    if target is None:
        target = "http://k-t.ep.interactivetool.localhost:8090/"

    def make_get_request(url):
        if "/entry_points?job_id=" in url:
            return _FakeResponse(entry_points)
        if "/jobs/JOB123?full=true" in url:
            return _FakeResponse({"state": job_state_value})
        if "/entry_points/ep1/access" in url:
            return _FakeResponse({"target": target})
        raise AssertionError(f"Unexpected GET {url}")

    def _post(payload):
        if launch_error is not None:
            raise launch_error
        return {"jobs": [{"id": "JOB123"}]}

    gi = mock.Mock()
    gi.url = "http://localhost:9090/api"
    gi.histories.create_history.return_value = {"id": "HIST1"}
    gi.tools._post.side_effect = _post
    gi.make_get_request.side_effect = make_get_request

    config = mock.Mock()
    config.user_gi = gi
    return config


def _serve(content="moo cow\n"):
    response = mock.Mock(text=content)
    response.raise_for_status = lambda: None
    return response


def test_run_interactive_tool_test_success_no_content():
    config = _make_config()
    with mock.patch("planemo.galaxy.test.interactive.requests.get", return_value=_serve()):
        result = run_interactive_tool_test(config, "interactivetool_simple")
    assert result["id"] == "interactivetool_simple-interactive"
    assert result["test_type"] == "interactive"
    assert result["has_data"] is True
    assert result["data"]["status"] == "success"
    assert result["data"]["content"] == "moo cow\n"
    assert result["data"]["output_problems"] == []
    assert result["data"]["execution_problem"] is None
    assert result["data"]["entry_points"] == [{"id": "ep1", "name": "Simple IT", "active": True}]
    assert result["data"]["job"]["id"] == "JOB123"


def test_run_interactive_tool_test_success_content_match():
    config = _make_config()
    with mock.patch("planemo.galaxy.test.interactive.requests.get", return_value=_serve()):
        result = run_interactive_tool_test(config, "interactivetool_simple", expected_content="moo cow\n")
    assert result["data"]["status"] == "success"


def test_run_interactive_tool_test_failure_content_mismatch():
    config = _make_config()
    with mock.patch("planemo.galaxy.test.interactive.requests.get", return_value=_serve("something else\n")):
        result = run_interactive_tool_test(config, "interactivetool_simple", expected_content="moo cow\n")
    assert result["data"]["status"] == "failure"
    assert any("does not match" in problem for problem in result["data"]["output_problems"])


def test_run_interactive_tool_test_error_launch_fails():
    config = _make_config(launch_error=RuntimeError("boom"))
    result = run_interactive_tool_test(config, "interactivetool_simple")
    assert result["data"]["status"] == "error"
    assert "Failed to launch" in result["data"]["execution_problem"]
    assert result["data"]["job"] is None


def test_run_interactive_tool_test_error_job_errors():
    config = _make_config(entry_points=[{"id": "ep1", "name": "Simple IT", "active": False}], job_state_value="error")
    result = run_interactive_tool_test(config, "interactivetool_simple", entry_point_timeout=5)
    assert result["data"]["status"] == "error"
    assert "did not become ready" in result["data"]["execution_problem"]


def test_run_interactive_tool_test_error_entry_points_time_out():
    config = _make_config(entry_points=[{"id": "ep1", "name": "Simple IT", "active": False}], job_state_value="running")
    with mock.patch("planemo.galaxy.test.interactive.time.sleep"):
        result = run_interactive_tool_test(config, "interactivetool_simple", entry_point_timeout=1)
    assert result["data"]["status"] == "error"
    assert "did not become ready" in result["data"]["execution_problem"]


def test_run_interactive_tool_test_error_content_unreachable():
    config = _make_config()
    with (
        mock.patch("planemo.galaxy.test.interactive.requests.get", side_effect=Exception("connection refused")),
        mock.patch("planemo.galaxy.test.interactive.time.sleep"),
    ):
        result = run_interactive_tool_test(config, "interactivetool_simple", content_timeout=1)
    assert result["data"]["status"] == "error"
    assert "did not respond" in result["data"]["execution_problem"]


def test_result_validates_as_report():
    config = _make_config()
    with mock.patch("planemo.galaxy.test.interactive.requests.get", return_value=_serve()):
        result = run_interactive_tool_test(config, "interactivetool_simple", expected_content="moo cow\n")
    report = PlanemoTestReport.model_validate({"version": "0.1", "tests": [result]})
    assert report.tests[0].data.status == "success"
    assert report.tests[0].test_type == "interactive"


def test_inputs_forwarded_to_tool_run():
    config = _make_config()
    with mock.patch("planemo.galaxy.test.interactive.requests.get", return_value=_serve()):
        run_interactive_tool_test(config, "interactivetool_simple", inputs={"param": "value"})
    posted = config.user_gi.tools._post.call_args.kwargs["payload"]
    assert posted["inputs"] == {"param": "value"}
    assert posted["tool_id"] == "interactivetool_simple"
