"""Unit tests for runnable test case checking and related functionality."""
import os

from planemo.runnable import (
    cases,
    for_path,
    SuccessfulRunResponse,
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


class MockRunResponse(SuccessfulRunResponse):
    start_datetime = None
    end_datetime = None

    def __init__(self, outputs_dict):
        self._outputs_dict = outputs_dict

    @property
    def log(self):
        return "My Log"

    @property
    def job_info(self):
        return {"command_line": "cat /tmp/1.txt > /tmp/2.out", "stdout": "Cat output", "stderr": "Cat problems..."}

    @property
    def outputs_dict(self):
        return self._outputs_dict

    @property
    def invocation_details(self):
        return None


__all__ = (
    "test_non_file_case_checker",
    "test_file_case_checker",
)
