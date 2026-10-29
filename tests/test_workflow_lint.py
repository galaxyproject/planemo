import json
import os
from importlib.metadata import version

import pytest
import yaml
from gxformat2.linting import Linter
from packaging.version import Version

from planemo.workflow_lint import (
    _lint_best_practices,
    _lint_schema_validation,
    _lint_structure,
    input_labels,
    WorkflowLintContext,
)
from .test_utils import TEST_DATA_DIR


def test_input_labels():
    wf_example_1 = os.path.join(TEST_DATA_DIR, "wf_repos", "basic_format2_ok", "basic_format2.gxwf.yml")
    labels = input_labels(wf_example_1)
    assert "the_input" in labels
    # verify non-input labels excluded
    assert "cat" not in labels

    wf_example_1 = os.path.join(TEST_DATA_DIR, "wf_repos", "basic_native_ok", "basic_native.yml")
    labels = input_labels(wf_example_1)
    assert "the_input" in labels
    # verify non-input labels excluded
    assert "cat" not in labels


@pytest.mark.parametrize("level", ["all", "error", "silent"])
def test_nested_native_diagnostics(tmp_path, capsys, level):
    workflow = {
        "a_galaxy_workflow": "true",
        "format-version": "0.1",
        "steps": {
            "0": {
                "type": "subworkflow",
                "subworkflow": {
                    "a_galaxy_workflow": "true",
                    "format-version": "0.1",
                    "steps": {"bad/~": {"type": "data_input", "label": "input"}},
                },
            }
        },
    }
    path = tmp_path / "nested.ga"
    path.write_text(json.dumps(workflow))
    context = WorkflowLintContext(level)
    context.lint("lint_structure", _lint_structure, str(path))

    assert context.failed("error")
    assert len(context.error_messages) == 1
    error = context.error_messages[0]
    assert error.message == "expected step_key to be integer not [bad/~]"
    assert error.linter == "NativeStepKeyNotInteger"
    assert len(context.warn_messages) == 2
    assert all(m.message == "Workflow contained no outputs" for m in context.warn_messages)
    # Older gxformat2 releases emitted these at the workflow root.
    if Version(version("gxformat2")) >= Version("0.29.0"):
        assert error.json_pointer == "/steps/0/subworkflow/steps/bad~1~0"
        assert [m.json_pointer for m in context.warn_messages] == ["/steps/0/subworkflow", ""]
    output = capsys.readouterr().out
    assert ("ERROR (NativeStepKeyNotInteger)" in output) == (level != "silent")
    assert ("WARNING" in output) == (level == "all")
    if level != "silent" and Version(version("gxformat2")) >= Version("0.29.0"):
        assert "[/steps/0/subworkflow/steps/bad~1~0]" in output


@pytest.mark.parametrize("suffix", [".ga", ".gxwf.yml"])
def test_workflow_training_topic(tmp_path, suffix):
    if suffix == ".ga":
        workflow = {
            "a_galaxy_workflow": "true",
            "format-version": "0.1",
            "steps": {"0": {"type": "data_input", "label": "input"}},
        }
    else:
        workflow = {"class": "GalaxyWorkflow", "steps": {"cat": {"tool_id": "cat1"}}}
    path = tmp_path / f"training{suffix}"
    path.write_text(yaml.safe_dump(workflow))
    context = WorkflowLintContext("silent")
    context.training_topic = "transcriptomics"
    context.lint("lint_structure", _lint_structure, str(path))
    assert "Missing tag(s)." in context.warn_messages
    assert "Missing workflow documentation (annotation or doc element)" in context.warn_messages
    assert context.failed("warn")
    assert not context.failed("error")


def test_format2_lint_phases(tmp_path, capsys):
    path = tmp_path / "invalid.gxwf.yml"
    path.write_text(
        yaml.safe_dump(
            {
                "class": "GalaxyWorkflow",
                "unexpected": True,
                "inputs": {},
                "outputs": {"out": {"outputSource": "absent/out"}},
                "steps": {"cat": {"tool_id": "cat1"}},
            }
        )
    )
    context = WorkflowLintContext("all")
    context.lint("lint_structure", _lint_structure, str(path))
    assert len(context.error_messages) == 1
    assert "references step 'absent'" in context.error_messages[0].message
    assert context.failed("error")
    capsys.readouterr()
    context.lint("lint_schema_validation", _lint_schema_validation, str(path))
    assert "unexpected" in capsys.readouterr().out
    context.lint("lint_best_practices", _lint_best_practices, str(path))
    output = capsys.readouterr().out
    assert "Workflow is not annotated." in output
    assert "references step 'absent'" not in output
    assert len(context.error_messages) == 1
    assert context.warn_messages


def test_gxformat2_context_aggregation(capsys):
    context = WorkflowLintContext("error", skip_types=["skipped"])
    context.training_topic = "transcriptomics"

    def lint_nested(target, lint_context):
        with lint_context.gxformat2_context() as workflow_context:
            child = workflow_context.child("steps").child("a/~").child("run")
            assert child.training_topic == "transcriptomics"
            child.warn("Nested warning", linter=Linter)
            child.error("Nested error", linter="NestedRule", json_pointer="/explicit")

    context.lint("lint_skipped", lint_nested, None)
    assert not context.message_list
    context.lint("lint_nested", lint_nested, None)
    assert context.failed("error")
    assert context.failed("warn")
    assert context.warn_messages[0].json_pointer == "/steps/a~1~0/run"
    assert context.warn_messages[0].linter == "Linter"
    assert context.error_messages[0].json_pointer == "/explicit"
    assert capsys.readouterr().out == ".. ERROR (NestedRule): Nested error [/explicit]\n"
    context.lint("lint_next", lambda target, ctx: ctx.info("Next phase"), None)
    assert len(context.message_list) == 3
    assert not capsys.readouterr().out
