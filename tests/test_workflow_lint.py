import os

from gxformat2.linting import Linter

from planemo.workflow_lint import (
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


def test_gxformat2_context_aggregation(capsys):
    context = WorkflowLintContext("error")
    context.training_topic = "transcriptomics"

    def lint_nested(target, lint_context):
        with lint_context.gxformat2_context() as workflow_context:
            child = workflow_context.child("steps").child("a/~").child("run")
            assert child.training_topic == "transcriptomics"
            child.warn("Nested warning", linter=Linter)
            child.error("Nested error", linter="NestedRule", json_pointer="/explicit")

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
