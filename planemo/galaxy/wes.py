"""A thin ``requests`` client for Galaxy's GA4GH WES endpoints.

This is a dependency-light wrapper (``requests`` only) covering the subset of the
WES wire protocol planemo needs to submit a workflow run and poll it to a
terminal state. It is modeled on the standalone ``gxy-wes`` reference client.
WES has no data-staging or output-download endpoint, so input staging and output
downloads are handled separately through the native Galaxy API.
"""

import json
from typing import (
    Any,
    Dict,
    Optional,
)

import requests

WES_PREFIX = "ga4gh/wes/v1"

STATE_COMPLETE = "COMPLETE"
# Terminal states that indicate the run did not succeed.
FAILURE_STATES = frozenset({"EXECUTOR_ERROR", "SYSTEM_ERROR", "CANCELED"})
# WES run states that mean "stop polling".
TERMINAL_STATES = FAILURE_STATES | {STATE_COMPLETE}


def is_terminal(state: str) -> bool:
    """Return True if ``state`` is a terminal WES run state (success or failure)."""
    return state in TERMINAL_STATES


def is_success(state: str) -> bool:
    """Return True if ``state`` is the successful terminal WES run state."""
    return state == STATE_COMPLETE


def is_failure(state: str) -> bool:
    """Return True if ``state`` is a terminal WES run state indicating failure."""
    return state in FAILURE_STATES


class WesError(Exception):
    """Raised when the WES server returns a non-2xx response."""

    def __init__(self, message: str, status_code: Optional[int] = None) -> None:
        super().__init__(message)
        #: HTTP status code of the failed response, when available.
        self.status_code = status_code


def detect_workflow_type(workflow_text: str) -> str:
    """Guess the Galaxy WES ``workflow_type`` from a workflow document.

    Returns ``gx_workflow_format2`` for Format2 (``class: GalaxyWorkflow``) or
    ``gx_workflow_ga`` for native ``.ga`` workflows. Galaxy re-validates the type
    against the referenced workflow, so this only needs to be approximately right.
    """
    text = workflow_text.lstrip()
    if text.startswith("{"):
        try:
            parsed = json.loads(workflow_text)
        except ValueError:
            parsed = {}
        if isinstance(parsed, dict) and parsed.get("class") == "GalaxyWorkflow":
            return "gx_workflow_format2"
        return "gx_workflow_ga"
    # YAML-ish: Format2 documents declare ``class: GalaxyWorkflow``.
    if "class: GalaxyWorkflow" in workflow_text or "class: 'GalaxyWorkflow'" in workflow_text:
        return "gx_workflow_format2"
    return "gx_workflow_ga"


class WesClient:
    """Minimal Galaxy GA4GH WES client."""

    def __init__(self, galaxy_url: str, api_key: Optional[str] = None, timeout: float = 60.0) -> None:
        self.galaxy_url = galaxy_url.rstrip("/")
        self.api_key = api_key
        self.timeout = timeout
        self.session = requests.Session()

    def _headers(self) -> Dict[str, str]:
        headers: Dict[str, str] = {}
        if self.api_key:
            headers["x-api-key"] = self.api_key
        return headers

    def _request(self, method: str, path: str, **kwargs: Any) -> requests.Response:
        url = f"{self.galaxy_url}/{path.lstrip('/')}"
        response = self.session.request(method, url, headers=self._headers(), timeout=self.timeout, **kwargs)
        if not response.ok:
            message = response.text
            try:
                body = response.json()
                message = body.get("err_msg", message)
            except ValueError:
                pass
            raise WesError(
                f"{method} {url} -> HTTP {response.status_code}: {message}", status_code=response.status_code
            )
        return response

    def submit_run(
        self,
        *,
        workflow_type: str,
        workflow_url: str,
        workflow_type_version: str = "1.0.0",
        params: Optional[Dict[str, Any]] = None,
        engine_parameters: Optional[Dict[str, Any]] = None,
    ) -> Dict[str, Any]:
        """Submit a workflow run, returning ``{"run_id": ...}``."""
        if not workflow_url:
            raise ValueError("workflow_url is required to submit a WES run")

        data: Dict[str, str] = {
            "workflow_type": workflow_type,
            "workflow_type_version": workflow_type_version,
            "workflow_url": workflow_url,
        }
        if params is not None:
            data["workflow_params"] = json.dumps(params)
        if engine_parameters is not None:
            data["workflow_engine_parameters"] = json.dumps(engine_parameters)
        return self._request("POST", f"{WES_PREFIX}/runs", data=data).json()

    def get_run_status(self, run_id: str) -> Dict[str, Any]:
        return self._request("GET", f"{WES_PREFIX}/runs/{run_id}/status").json()


__all__ = (
    "detect_workflow_type",
    "FAILURE_STATES",
    "is_failure",
    "is_success",
    "is_terminal",
    "STATE_COMPLETE",
    "TERMINAL_STATES",
    "WesClient",
    "WesError",
)
