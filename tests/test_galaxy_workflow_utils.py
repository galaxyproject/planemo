"""Test utilities for dealing with Galaxy workflows."""

import os

from planemo.galaxy.workflows import (
    describe_outputs,
    output_stubs_for_workflow,
    required_input_labels,
)
from planemo.runnable import for_path
from .test_utils import TEST_DATA_DIR


def test_describe_outputs():
    wf_path = os.path.join(TEST_DATA_DIR, "wf1.gxwf.yml")
    runnable = for_path(wf_path)
    outputs = describe_outputs(runnable)
    assert len(outputs) == 1
    output = outputs[0]
    assert output.order_index == 1
    assert output.output_name == "out_file1"
    assert output.label == "wf_output_1"


def test_required_input_steps_excludes_falsy_defaults():
    """Inputs with a default are not required, even when the default is falsy (boolean false, int 0)."""
    wf_path = os.path.join(TEST_DATA_DIR, "wf_required_defaults.gxwf.yml")
    required = set(required_input_labels(wf_path))
    assert "needs_value" in required
    assert "bool_default_false" not in required
    assert "int_default_zero" not in required


def test_output_stubs_skip_anonymous_and_empty_labels():
    """Empty/anonymous output labels are filtered via gxformat2's canonical predicate."""
    wf_path = os.path.join(TEST_DATA_DIR, "wf_empty_output_label.gxwf.yml")
    stubs = output_stubs_for_workflow(wf_path)
    assert stubs == {}
