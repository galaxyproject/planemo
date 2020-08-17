"""Unit tests for runnable test case checking and related functionality."""
import os

from planemo.runnable import (
    cases,
    for_path,
)
from .test_utils import TEST_DATA_DIR


def test_non_file_case_checker():
    """Verify simply usage of :func:`planemo.runnable.TestCase.check`."""
    int_tool_path = os.path.join(TEST_DATA_DIR, "int_tool.cwl")
    test_cases = cases(for_path(int_tool_path))
    assert len(test_cases) == 1
    test_case = test_cases[0]
    outputs_dict = {
        "output": 4,
    }

    sd = test_case.structured_test_data(MockRunResponse(outputs_dict))
    assert sd["data"]["status"] == "success"

    bad_outputs_dict = {
        "output": 5,
    }
    sd = test_case.structured_test_data(MockRunResponse(bad_outputs_dict))
    assert sd["data"]["status"] == "failure"


def test_file_case_checker():
    hello_txt_path = os.path.join(TEST_DATA_DIR, "hello.txt")
    int_tool_path = os.path.join(TEST_DATA_DIR, "cat_tool.cwl")
    test_cases = cases(for_path(int_tool_path))
    assert len(test_cases) == 1
    test_case = test_cases[0]
    outputs_dict = {
        "output_file": {
            "path": hello_txt_path,
        }
    }

    sd = test_case.structured_test_data(MockRunResponse(outputs_dict))
    assert sd["data"]["status"] == "success"

    not_hello_txt_path = os.path.join(TEST_DATA_DIR, "int_tool_job.json")
    bad_outputs_dict = {
        "output_file": {
            "path": not_hello_txt_path,
        }
    }
    sd = test_case.structured_test_data(MockRunResponse(bad_outputs_dict))
    assert sd["data"]["status"] == "failure"


class MockRunResponse(object):

    def __init__(self, outputs_dict):
        self.outputs_dict = outputs_dict
        self.was_successful = True
        self.job_info = {
            "command_line": "cat /tmp/1.txt > /tmp/2.out",
            "stdout": "Cat output",
            "stderr": "Cat problems..."
        }
        self.invocation_details = None
        self.log = "My Log"


class MockErrorRunResponse(object):
    was_successful = False
    error_message = "Problem running job."
    log = "Error log"


__all__ = (
    "test_non_file_case_checker",
    "test_file_case_checker",
)
