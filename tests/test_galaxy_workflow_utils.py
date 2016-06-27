"""Test utilities for dealing with Galaxy workflows."""
import os

from .test_utils import TEST_DATA_DIR
from planemo.galaxy.workflows import describe_outputs


def test_describe_outputs():
    wf_path = os.path.join(TEST_DATA_DIR, "wf1.gxwf.yml")
    outputs = describe_outputs(wf_path)
    assert len(outputs) == 1
    output = outputs[0]
    assert output.order_index == 1
    assert output.output_name == "out_file1"
    assert output.label == "wf_output_1"
