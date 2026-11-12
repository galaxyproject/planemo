"""Actions for testing Galaxy Interactive Tools (GxITs).

An interactive tool, when run, launches a container that opens one or more
network entry points (ports). "Testing" such a tool means verifying that it
starts up and that its entry points become active and reachable through the
Galaxy Interactive Tools (GxIT) proxy.

The logic here mirrors the integration test added in planemo PR #1568
("Enable GxITs in planemo serve") so that the same behaviour is available as a
first class ``planemo`` feature.
"""

import time
from datetime import datetime
from typing import (
    Any,
    Dict,
    List,
    Optional,
    Tuple,
)

import requests
from galaxy.tool_util.parser import get_tool_source


def is_interactive_tool(path: str) -> bool:
    """Return True iff the tool defined at ``path`` is an interactive tool."""
    try:
        tool_source = get_tool_source(path)
    except Exception:
        return False
    return tool_source.parse_tool_type() == "interactive"


def tool_id_for_path(path: str) -> str:
    """Return the tool id for the tool defined at ``path``."""
    return get_tool_source(path).parse_id()


def test_case_count(path: str) -> int:
    """Return the number of embedded ``<tests>`` cases declared by the tool at ``path``."""
    return len(get_tool_source(path).parse_tests_to_dict().get("tests", []))


def inputs_from_test_case(path: str, test_index: int) -> Dict[str, Any]:
    """Extract the tool inputs declared in ``test_index`` (1-indexed) of the tool at ``path``.

    Interactive tools that require parameters can be tested by reusing the
    inputs of one of the tool's embedded ``<tests>`` cases, mirroring how
    ``planemo test`` selects a test with ``--test_index`` (counting starts at 1).

    :raises ValueError: if ``test_index`` is not a valid 1-based index into the
        tool's embedded test cases.
    """
    tool_source = get_tool_source(path)
    tests = tool_source.parse_tests_to_dict().get("tests", [])
    
    if not 1 <= test_index <= len(tests):
        raise ValueError(
            f"Tool [{tool_source.parse_id()}] has {len(tests)} test case(s); "
            f"test_index [{test_index}] is out of range."
        )
    test_case = tests[test_index - 1]
    return {inp["name"]: inp["value"] for inp in test_case.get("inputs", [])}


def resolve_interactive_inputs(
    path: str,
    test_index: Optional[int],
    manual_inputs: Optional[Dict[str, Any]],
) -> Dict[str, Any]:
    """Resolve the inputs to launch the tool at ``path`` with.

    The inputs of the referenced embedded test case (``test_index``) form the
    base; any explicit ``manual_inputs`` override them. Returns an empty dict
    when neither is supplied.
    """
    inputs: Dict[str, Any] = {}
    if test_index is not None:
        inputs.update(inputs_from_test_case(path, test_index))
    if manual_inputs:
        inputs.update(manual_inputs)
    return inputs


def launch_interactive_tool(user_gi, tool_id: str, inputs: Optional[Dict[str, Any]] = None) -> str:
    """Create a history and run ``tool_id`` in it, returning the job id.

    :param user_gi: a bioblend :class:`~bioblend.galaxy.GalaxyInstance`
    :param tool_id: the interactive tool id to run
    :param inputs: tool inputs (default: empty)
    :return: the id of the job created by the tool run
    """
    if inputs is None:
        inputs = {}
    history = user_gi.histories.create_history("Interactive Tool Test")
    tool_run = user_gi.tools._post(
        payload={
            "tool_id": tool_id,
            "history_id": history["id"],
            "inputs": inputs,
        }
    )
    if "jobs" not in tool_run:
        raise RuntimeError(f"Running interactive tool [{tool_id}] did not return a job: {tool_run}")
    return tool_run["jobs"][0]["id"]


def job_state(gi, job_id: str) -> Optional[str]:
    """Return the state of ``job_id`` (best effort, ``None`` on failure)."""
    try:
        resp = gi.make_get_request(f"{gi.url}/jobs/{job_id}?full=true")
        resp.raise_for_status()
        return resp.json().get("state")
    except Exception:
        return None


def wait_for_active_entry_points(gi, job_id: str, timeout: int = 120) -> List[Dict[str, Any]]:
    """Poll ``GET /api/entry_points?job_id=...`` until all entry points are active.

    Raises a :class:`TimeoutError` if the entry points do not become active
    within ``timeout`` seconds and an :class:`Exception` if the job errors.
    """
    end = time.time() + timeout
    while time.time() < end:
        resp = gi.make_get_request(f"{gi.url}/entry_points?job_id={job_id}")
        resp.raise_for_status()
        entry_points = resp.json()
        if entry_points and all(ep.get("active") for ep in entry_points):
            return entry_points
        # Check if job errored
        state = job_state(gi, job_id)
        if state == "error":
            raise RuntimeError(f"Interactive tool job {job_id} failed (state: error).")
        time.sleep(3)
    raise TimeoutError(f"Entry points for job {job_id} did not become active within {timeout}s.")


def get_entry_point_target(gi, entry_point_id: str) -> str:
    """Return the GxIT proxy target URL for an entry point."""
    resp = gi.make_get_request(f"{gi.url}/entry_points/{entry_point_id}/access")
    resp.raise_for_status()
    return resp.json()["target"]


def _split_proxied_url(target_url: str) -> Tuple[str, Optional[str]]:
    """Split a GxIT proxy target URL into ``(request_url, host_header)``.

    The target URL has the form
    ``http://{key}-{token}.ep.interactivetool.{host}:{port}{path}``. The GxIT
    proxy is addressed by its host/port with the entry-point subdomain supplied
    via the ``Host`` header, so the subdomain is parsed out and returned as the
    ``Host`` header value. When no ``.interactivetool.`` marker is present the
    URL is returned unchanged with no ``Host`` header.
    """
    scheme, rest = target_url.split("://", 1)
    host_header = rest.split("/", 1)[0] if "/" in rest else rest
    if ".interactivetool." in rest:
        _prefix, host_and_port = rest.split(".interactivetool.", 1)
        # host_and_port is like "localhost:8090/" - strip any trailing path
        host_and_port = host_and_port.split("/", 1)[0]
        return f"{scheme}://{host_and_port}", host_header
    return target_url, None


def fetch_proxied_content(target_url: str, timeout: int = 30) -> str:
    """Request content through the GxIT proxy for an entry point target URL.

    Retries until ``timeout`` seconds elapse; raises :class:`TimeoutError` on
    failure. See :func:`_split_proxied_url` for how the target URL is resolved.
    """
    end = time.time() + timeout
    last_error: Optional[Exception] = None
    while time.time() < end:
        try:
            url, host_header = _split_proxied_url(target_url)
            headers = {"Host": host_header} if host_header else {}
            response = requests.get(url, timeout=5, headers=headers)
            response.raise_for_status()
            return response.text
        except Exception as e:
            last_error = e
            time.sleep(2)
    raise TimeoutError(f"Failed to get proxied content from {target_url} within {timeout}s: {last_error}")


def run_interactive_tool_test(
    config,
    tool_id: str,
    inputs: Optional[Dict[str, Any]] = None,
    expected_content: Optional[str] = None,
    entry_point_timeout: int = 120,
    content_timeout: int = 30,
) -> Dict[str, Any]:
    """Run the interactive-tool smoke test and return a structured test case.

    The test launches the tool, waits for its entry points to become active,
    fetches content through the GxIT proxy and (optionally) asserts it matches
    ``expected_content``.

    :param config: a planemo Galaxy config exposing a ``user_gi`` attribute
    :param tool_id: the interactive tool id to test
    :param inputs: tool inputs (default: empty)
    :param expected_content: when set, the served content of the primary
        (first) entry point must match this exactly
    :param entry_point_timeout: seconds to wait for entry points to activate
    :param content_timeout: seconds to wait for proxied content
    :return: a planemo structured test case dict (see ``PlanemoTestCase``)
    """
    if inputs is None:
        inputs = {}
    start = datetime.now()
    data: Dict[str, Any] = {
        "status": "error",
        "inputs": inputs,
        "job": None,
        "entry_points": None,
        "content": None,
        "output_problems": [],
        "execution_problem": None,
        "problem_log": None,
        "start_datetime": start.isoformat(),
        "end_datetime": None,
    }
    test_case: Dict[str, Any] = {
        "id": f"{tool_id}-interactive",
        "has_data": True,
        "test_type": "interactive",
        "data": data,
    }

    def finish(status: str) -> Dict[str, Any]:
        data["status"] = status
        data["end_datetime"] = datetime.now().isoformat()
        return test_case

    user_gi = config.user_gi

    # 1. Launch the tool.
    try:
        job_id = launch_interactive_tool(user_gi, tool_id, inputs)
    except Exception as e:
        data["execution_problem"] = f"Failed to launch interactive tool [{tool_id}]: {e}"
        return finish("error")
    data["job"] = {"id": job_id, "state": None}

    # 2. Wait for the entry points to become active.
    try:
        entry_points = wait_for_active_entry_points(user_gi, job_id, timeout=entry_point_timeout)
    except Exception as e:
        data["execution_problem"] = f"Interactive tool [{tool_id}] did not become ready: {e}"
        data["job"]["state"] = job_state(user_gi, job_id)
        return finish("error")
    data["entry_points"] = entry_points
    data["job"]["state"] = job_state(user_gi, job_id)

    # 3. Fetch content through the proxy to confirm the port is actually serving.
    entry_point = entry_points[0]
    entry_point_name = entry_point.get("name") or entry_point.get("id")
    try:
        target = get_entry_point_target(user_gi, entry_point["id"])
        content = fetch_proxied_content(target, timeout=content_timeout)
    except Exception as e:
        data["execution_problem"] = f"Entry point [{entry_point_name}] did not respond: {e}"
        return finish("error")
    data["content"] = content

    # 4. Optionally assert the served content matches the expected content.
    if expected_content is not None and content != expected_content:
        data["status"] = "failure"
        data["output_problems"].append(
            f"Entry point [{entry_point_name}] content [{content!r}] does not match "
            f"expected content [{expected_content!r}]."
        )
        return finish("failure")

    return finish("success")


__all__ = (
    "fetch_proxied_content",
    "get_entry_point_target",
    "inputs_from_test_case",
    "is_interactive_tool",
    "job_state",
    "launch_interactive_tool",
    "run_interactive_tool_test",
    "resolve_interactive_inputs",
    "test_case_count",
    "tool_id_for_path",
    "wait_for_active_entry_points",
)
