"""Test utilities for dealing with Galaxy workflows."""
import os

from planemo.galaxy.workflows import describe_outputs
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
