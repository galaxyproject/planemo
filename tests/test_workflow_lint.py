import os

from planemo.workflow_lint import input_labels
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
