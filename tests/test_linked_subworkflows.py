"""Tests for workflows whose subworkflow steps are linked (URL / TRS references) rather than embedded."""

import base64
import json
import os

import yaml

from planemo import autoupdate
from planemo.galaxy import workflows
from planemo.galaxy.workflows import (
    expand_linked_subworkflows,
    export_workflow_dict,
    find_tool_ids,
    get_tool_ids_for_workflow,
    linked_subworkflow_url,
    load_shed_repos,
)
from planemo.runnable import for_path
from planemo.workflow_lint import (
    _iter_tool_steps,
    _lint_linked_subworkflows,
    _load_workflow_dict,
    WorkflowLintContext,
)
from .test_utils import TEST_DATA_DIR

PARENT_TOOL_ID = "toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_head_tool/9.3+galaxy1"
CHILD_TOOL_ID = "toolshed.g2.bx.psu.edu/repos/devteam/fastqc/fastqc/0.72+galaxy1"
TRS_URL = (
    "https://dockstore.org/api/ga4gh/trs/v2/tools/#workflow/github.com/iwc-workflows/linked-child/main/versions/v0.1"
)
TRS_DESCRIPTOR_URL = (
    "https://dockstore.org/api/ga4gh/trs/v2/tools/"
    "%23workflow%2Fgithub.com%2Fiwc-workflows%2Flinked-child%2Fmain/versions/v0.1/GALAXY/descriptor"
)
URL_GA = os.path.join(TEST_DATA_DIR, "wf_linked_subworkflow_url.ga")
TRS_GA = os.path.join(TEST_DATA_DIR, "wf_linked_subworkflow_trs.ga")
URL_GXWF = os.path.join(TEST_DATA_DIR, "wf_linked_subworkflow_url.gxwf.yml")


def _child_content():
    with open(os.path.join(TEST_DATA_DIR, "wf_linked_subworkflow_child.ga")) as f:
        return f.read()


class _FakeResponse:
    def __init__(self, payload):
        self._payload = payload

    def raise_for_status(self):
        pass

    def json(self):
        return self._payload


def _fake_trs(monkeypatch):
    """Serve the child workflow from a fake Dockstore TRS descriptor endpoint."""
    requested = []

    def fake_get(url, timeout=None):
        requested.append(url)
        assert url == TRS_DESCRIPTOR_URL
        return _FakeResponse({"content": _child_content()})

    monkeypatch.setattr(workflows.requests, "get", fake_get)
    return requested


def _load(path):
    with open(path) as f:
        return yaml.safe_load(f)


def test_linked_subworkflow_url():
    b64 = "base64://" + base64.b64encode(b"{}").decode()
    assert linked_subworkflow_url({"type": "subworkflow", "content_source": "url", "content_id": b64}) == b64
    assert (
        linked_subworkflow_url({"type": "subworkflow", "content_source": "trs_url", "content_id": TRS_URL}) == TRS_URL
    )
    # gxformat2-converted steps carry the URL in content_id without content_source
    assert linked_subworkflow_url({"type": "subworkflow", "content_id": TRS_URL}) == TRS_URL
    trs_ids_step = {
        "type": "subworkflow",
        "trs_server": "dockstore",
        "trs_tool_id": "#workflow/github.com/iwc-workflows/linked-child/main",
        "trs_version_id": "v0.1",
    }
    assert linked_subworkflow_url(trs_ids_step) == TRS_URL
    assert linked_subworkflow_url({"run": TRS_URL}) == TRS_URL
    # embedded subworkflows, local imports and tool steps aren't linked
    assert linked_subworkflow_url({"type": "subworkflow", "subworkflow": {"steps": {}}}) is None
    assert linked_subworkflow_url({"run": "child.gxwf.yml"}) is None
    assert linked_subworkflow_url({"type": "tool", "content_id": CHILD_TOOL_ID}) is None


def test_get_tool_ids_for_workflow_skips_linked_subworkflows():
    # autoupdate must only consider the workflow's own tools, linked subworkflows are pinned references
    for path in (URL_GA, TRS_GA, URL_GXWF):
        assert get_tool_ids_for_workflow(_load(path)) == [PARENT_TOOL_ID]


def test_find_tool_ids_expands_url_linked_subworkflows():
    for path in (URL_GA, URL_GXWF):
        assert find_tool_ids(path) == [PARENT_TOOL_ID, CHILD_TOOL_ID]


def test_find_tool_ids_expands_trs_linked_subworkflows(monkeypatch):
    requested = _fake_trs(monkeypatch)
    assert find_tool_ids(TRS_GA) == [PARENT_TOOL_ID, CHILD_TOOL_ID]
    assert requested == [TRS_DESCRIPTOR_URL]


def test_expand_trs_id_linked_subworkflow(monkeypatch):
    _fake_trs(monkeypatch)
    workflow = _load(TRS_GA)
    step = workflow["steps"]["2"]
    del step["content_source"], step["content_id"]
    step.update(
        trs_server="dockstore",
        trs_tool_id="#workflow/github.com/iwc-workflows/linked-child/main",
        trs_version_id="v0.1",
    )
    expanded = expand_linked_subworkflows(workflow)
    assert get_tool_ids_for_workflow(expanded) == [PARENT_TOOL_ID, CHILD_TOOL_ID]
    # the input is left untouched
    assert "subworkflow" not in workflow["steps"]["2"]


def test_expand_without_linked_subworkflows_is_noop():
    workflow = _load(os.path.join(TEST_DATA_DIR, "wf_same_name_outputs.ga"))
    assert expand_linked_subworkflows(workflow) is workflow


def test_embedded_gxformat2_subworkflow_with_explicit_type(tmp_path):
    # an embedded run with an explicit "type: subworkflow" is not a linked subworkflow
    workflow = {
        "class": "GalaxyWorkflow",
        "steps": {
            "nested": {
                "type": "subworkflow",
                "run": {"class": "GalaxyWorkflow", "steps": {"cat": {"tool_id": CHILD_TOOL_ID}}},
            }
        },
    }
    assert get_tool_ids_for_workflow(workflow) == [CHILD_TOOL_ID]
    assert [step["tool_id"] for step in _iter_tool_steps(workflow)] == [CHILD_TOOL_ID]
    path = tmp_path / "wf.gxwf.yml"
    path.write_text(yaml.safe_dump(workflow))
    assert _load_workflow_dict(str(path))["steps"]["nested"]["run"]["class"] == "GalaxyWorkflow"
    lint_context = WorkflowLintContext("all")
    _lint_linked_subworkflows(str(path), lint_context)
    assert not lint_context.found_errors


def test_load_shed_repos_includes_linked_subworkflow_tools():
    repos = load_shed_repos(for_path(URL_GA))
    assert [(r["name"], r["revisions"]) for r in repos] == [
        ("text_processing", ["ddf54b12c295"]),
        ("fastqc", ["e7b2202befea"]),
    ]
    repos = load_shed_repos(for_path(URL_GXWF))
    assert [(r["name"], r["changeset_revision"]) for r in repos] == [
        ("text_processing", "ddf54b12c295"),
        ("fastqc", "e7b2202befea"),
    ]
    assert all(r["tool_shed_url"] == "https://toolshed.g2.bx.psu.edu" for r in repos)


def test_load_shed_repos_trs_linked_subworkflow(monkeypatch):
    _fake_trs(monkeypatch)
    repos = load_shed_repos(for_path(TRS_GA))
    assert [r["name"] for r in repos] == ["text_processing", "fastqc"]


def test_lint_skips_linked_subworkflow_tool_steps():
    for path in (URL_GA, TRS_GA, URL_GXWF):
        assert [step["tool_id"] for step in _iter_tool_steps(_load(path))] == [PARENT_TOOL_ID]


def test_lint_workflow_dict_drops_linked_references():
    # keeps gxformat2 lint from fetching linked subworkflows or flagging content_source in its strict schema
    step = _load_workflow_dict(TRS_GA)["steps"]["2"]
    assert "content_source" not in step and "content_id" not in step
    assert "run" not in _load_workflow_dict(URL_GXWF)["steps"]["linked subworkflow"]


def _lint_linked(tmp_path, content_source, content_id):
    workflow = _load(TRS_GA)
    workflow["steps"]["2"].update(content_source=content_source, content_id=content_id)
    path = tmp_path / "wf.ga"
    path.write_text(json.dumps(workflow))
    lint_context = WorkflowLintContext("all")
    _lint_linked_subworkflows(str(path), lint_context)
    return lint_context


def test_lint_linked_subworkflows(tmp_path):
    lint_context = _lint_linked(tmp_path, "trs_url", TRS_URL)
    assert not lint_context.found_errors and not lint_context.found_warns
    assert "pinned to TRS version 'v0.1'" in lint_context.valid_messages[0].message

    lint_context = _lint_linked(tmp_path, "trs_url", TRS_URL.replace("versions/v0.1", "versions/main"))
    assert "looks like a branch" in lint_context.warn_messages[0].message

    lint_context = _lint_linked(tmp_path, "trs_url", TRS_URL.replace("/versions/v0.1", ""))
    assert "is not a valid TRS URL" in lint_context.error_messages[0].message

    lint_context = _lint_linked(tmp_path, "url", "https://example.org/child.ga")
    assert not lint_context.found_errors and not lint_context.found_warns
    assert "not pinned" in lint_context.info_messages[0].message

    lint_context = _lint_linked(tmp_path, None, None)
    assert "neither embeds a subworkflow nor references one" in lint_context.error_messages[0].message


class _FakeGalaxy:
    url = "https://galaxy.example/api"

    def __init__(self, workflows_by_id):
        self._workflows_by_id = workflows_by_id
        self.requests = []
        self.refactored = []
        self.workflows = self

    def make_get_request(self, url, params=None):
        self.requests.append((url, params))
        return _FakeResponse(self._workflows_by_id.get(url.rsplit("/", 1)[-1]))

    def refactor_workflow(self, workflow_id, actions):
        self.refactored.append(workflow_id)


def test_autoupdate_does_not_upgrade_linked_subworkflows():
    gi = _FakeGalaxy(
        {
            "parent": {
                "id": "parent",
                "steps": {
                    "0": {"type": "subworkflow", "workflow_id": "embedded"},
                    "1": {"type": "subworkflow", "workflow_id": "linked"},
                },
            },
            "embedded": {"id": "embedded", "steps": {}},
            "linked": {"id": "linked", "source_metadata": {"trs_url": TRS_URL}, "steps": {}},
        }
    )

    class FakeConfig:
        user_gi = gi

    autoupdate._update_wf(FakeConfig(), "parent")
    assert gi.refactored == ["embedded", "parent"]


def test_export_workflow_dict_preserves_linked_subworkflows():
    gi = _FakeGalaxy({"download": {"name": "exported"}})
    assert export_workflow_dict(gi, "wf_id") == {"name": "exported"}
    assert gi.requests == [
        ("https://galaxy.example/api/workflows/wf_id/download", {"preserve_external_subworkflow_links": "true"})
    ]
