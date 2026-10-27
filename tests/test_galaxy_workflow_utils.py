"""Test utilities for dealing with Galaxy workflows."""

import os

from planemo.galaxy.workflows import (
    describe_outputs,
    job_template,
    job_template_with_metadata,
    output_stubs_for_workflow,
    parameter_input_labels,
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


def test_describe_outputs_dict_tool_state():
    """A .ga with dict (not JSON-string) tool_state must not crash describe_outputs."""
    wf_path = os.path.join(TEST_DATA_DIR, "wf14-unlinted-best-practices-dict-tool-state.ga")
    runnable = for_path(wf_path)
    outputs = describe_outputs(runnable)
    assert len(outputs) == 1
    assert outputs[0].output_name == "outfile"


def test_required_input_steps_excludes_falsy_defaults():
    """Inputs with a default are not required, even when the default is falsy (boolean false, int 0)."""
    wf_path = os.path.join(TEST_DATA_DIR, "wf_required_defaults.gxwf.yml")
    required = set(required_input_labels(wf_path))
    assert "needs_value" in required
    assert "bool_default_false" not in required
    assert "int_default_zero" not in required


def test_parameter_input_labels():
    """Parameter inputs, including multiple-value ones, but not data or collection inputs."""
    wf_path = os.path.join(TEST_DATA_DIR, "wf_multiple_param_inputs.gxwf.yml")
    assert parameter_input_labels(workflow_path=wf_path) == ["text_list", "int_list"]
    wf_path = os.path.join(TEST_DATA_DIR, "wf_param_types.gxwf.yml")
    assert parameter_input_labels(workflow_path=wf_path) == [
        "int_param",
        "int_with_default",
        "string_param",
        "bool_param",
    ]


def test_parameter_input_labels_native():
    wf_path = os.path.join(TEST_DATA_DIR, "wf_native_param_inputs.ga")
    assert parameter_input_labels(workflow_path=wf_path) == ["text_list", "int_param"]


def test_parameter_input_labels_untyped_input_is_data():
    wf_path = os.path.join(TEST_DATA_DIR, "wf1.gxwf.yml")
    assert parameter_input_labels(workflow_path=wf_path) == []


def test_output_stubs_skip_anonymous_and_empty_labels():
    """Empty/anonymous output labels are filtered via gxformat2's canonical predicate."""
    wf_path = os.path.join(TEST_DATA_DIR, "wf_empty_output_label.gxwf.yml")
    stubs = output_stubs_for_workflow(wf_path)
    assert stubs == {}


FILE_TEMPLATE = {"class": "File", "path": "todo_test_data_path.ext"}
LIST_TEMPLATE = {
    "class": "Collection",
    "collection_type": "list",
    "elements": [{"class": "File", "identifier": "todo_element_name", "path": "todo_test_data_path.ext"}],
}


def test_job_template_input_kinds():
    """Data, collection and scalar parameter inputs get templates; multiple-value parameters a TODO."""
    wf_path = os.path.join(TEST_DATA_DIR, "wf_native_param_inputs.ga")
    expected = {
        "text_list": {"TODO"},
        "int_param": "todo_param_value",
        "data_input": FILE_TEMPLATE,
        "data_list": LIST_TEMPLATE,
    }
    assert job_template(wf_path) == expected
    template, metadata = job_template_with_metadata(wf_path)
    assert template == expected
    assert {label: m["type"] for label, m in metadata.items()} == {
        "text_list": ["string"],
        "int_param": "int",
        "data_input": "data",
        "data_list": "collection",
    }


def test_job_template_untyped_input_is_data():
    wf_path = os.path.join(TEST_DATA_DIR, "wf1.gxwf.yml")
    assert job_template(wf_path) == {"input1": FILE_TEMPLATE}
    template, _ = job_template_with_metadata(wf_path)
    assert template == {"input1": FILE_TEMPLATE}
