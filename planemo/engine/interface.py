"""Module contianing the :class:`Engine` abstraction."""

import abc
import copy
import json
import os
import tempfile
from contextlib import contextmanager
from typing import (
    Any,
    Callable,
    Dict,
    Iterator,
    List,
    Optional,
)
from urllib.parse import urlparse

import click

from planemo.exit_codes import EXIT_CODE_UNSUPPORTED_FILE_TYPE
from planemo.io import error
from planemo.runnable import (
    cases,
    RunnableType,
    TestCase,
)
from planemo.test.results import StructuredData


def _absolute_test_data_path(path: Any, tests_directory: str) -> Any:
    if not isinstance(path, str) or os.path.isabs(path) or urlparse(path).scheme or path.startswith("#"):
        return path
    return os.path.abspath(os.path.join(tests_directory, path))


def _absolutize_composite_data(file_value: Dict[str, Any], tests_directory: str) -> None:
    composite_data = file_value.get("composite_data") or []
    for composite_index, composite_item in enumerate(composite_data):
        if isinstance(composite_item, dict):
            for path_key in ("path", "location"):
                if path_key in composite_item:
                    composite_item[path_key] = _absolute_test_data_path(composite_item[path_key], tests_directory)
        elif isinstance(composite_item, str):
            composite_data[composite_index] = _absolute_test_data_path(composite_item, tests_directory)


def _absolutize_nested_job_paths(value: Any, tests_directory: str) -> None:
    if isinstance(value, list):
        for item in value:
            _absolutize_nested_job_paths(item, tests_directory)
    elif isinstance(value, dict):
        if value.get("class") in ("File", "Directory"):
            for path_key in ("path", "location"):
                if path_key in value:
                    value[path_key] = _absolute_test_data_path(value[path_key], tests_directory)
            _absolutize_composite_data(value, tests_directory)

        for item in value.values():
            _absolutize_nested_job_paths(item, tests_directory)


def _absolutize_job_paths(job: Dict[str, Any], tests_directory: str) -> Dict[str, Any]:
    """Copy a test job and resolve its local data paths against the test directory."""
    prepared_job = copy.deepcopy(job)
    _absolutize_nested_job_paths(prepared_job, tests_directory)
    return prepared_job


@contextmanager
def materialized_job_paths(test_cases: List[TestCase]) -> Iterator[List[str]]:
    """Yield a job file path for each test case, for the duration of the context.

    Jobs defined inline in a test definition get written to a temporary directory
    rather than beside the definition, so the source directory is left untouched and
    need not be writable. Their relative data paths are resolved against the test
    definition's directory first so they survive the move - the test case's own job
    is not modified. Test cases that already point at a job file are passed through
    untouched.
    """
    with tempfile.TemporaryDirectory(prefix="planemo-test-jobs-") as job_directory:
        job_paths = []
        for index, test_case in enumerate(test_cases):
            if test_case.job_path is not None:
                job_paths.append(test_case.job_path)
                continue
            # a test case defines exactly one of job_path and job
            assert test_case.job is not None
            job_path = os.path.join(job_directory, f"job-{index}.json")
            with open(job_path, "w") as f:
                json.dump(_absolutize_job_paths(test_case.job, test_case.tests_directory), f)
            job_paths.append(job_path)
        yield job_paths


class Engine(metaclass=abc.ABCMeta):
    """Abstract description of an external process for running tools or workflows."""

    @abc.abstractmethod
    def run(self, path, job_path):
        """Run a job using a compatible artifact (workflow or tool)."""

    @abc.abstractmethod
    def cleanup(self):
        """Release any resources used to run/test with this engine."""

    @abc.abstractmethod
    def test(self, runnables):
        """Test runnable artifacts (workflow or tool)."""


class BaseEngine(Engine):
    """Base class providing context and keywords for Engine implementations."""

    handled_runnable_types: List[RunnableType] = []

    def __init__(self, ctx, **kwds):
        """Store context and kwds."""
        self._ctx = ctx
        self._kwds = kwds

    def can_run(self, runnable):
        """Use subclass's ``handled_runnable_types`` variable to infer ``can_run``."""
        return runnable.type in self.handled_runnable_types

    def cleanup(self):
        """Default no-op cleanup method."""

    def run(self, runnables, job_paths, output_collectors: Optional[List[Callable]] = None):
        """Run a job using a compatible artifact (workflow or tool)."""
        self._check_can_run_all(runnables)
        run_responses = self._run(runnables, job_paths, output_collectors)
        return run_responses

    @abc.abstractmethod
    def _run(
        self,
        runnables,
        job_path,
        output_collectors: Optional[List[Callable]] = None,
        test_timeout: Optional[int] = None,
    ):
        """Run a job using a compatible artifact (workflow or tool) wrapped as a runnable."""

    def _check_can_run(self, runnable):
        if not self.can_run(runnable):
            template = "Engine type [%s] cannot execute [%s]s"
            message = template % (self.__class__, runnable.type)
            error(message)
            self._ctx.exit(EXIT_CODE_UNSUPPORTED_FILE_TYPE)

    def _check_can_run_all(self, runnables):
        for runnable in runnables:
            self._check_can_run(runnable)

    def test(self, runnables, test_timeout):
        """Test runnable artifacts (workflow or tool)."""
        self._check_can_run_all(runnables)
        test_cases = [t for tl in map(cases, runnables) for t in tl]

        # Filter test cases by specified indices if provided
        test_indices = self._kwds.get("test_index", ())
        if any(i < 1 for i in test_indices):
            raise ValueError("test_index must be 1-based (>= 1)")
        if test_indices:
            filtered_test_cases = [tc for i, tc in enumerate(test_cases, start=1) if i in test_indices]
            if filtered_test_cases:
                test_cases = filtered_test_cases
            else:
                # If no tests match the specified indices, log a warning and use original
                self._ctx.log(f"Warning: No tests found with indices {test_indices}. Running all tests instead.")

        # Filter test cases to only previously-failed ones when --failed/--lf is set
        if self._kwds.get("failed"):
            failed_json = self._kwds.get("failed_json") or self._kwds.get("test_output_json")
            if not failed_json or not os.path.exists(failed_json):
                raise click.ClickException(
                    "--failed/--lf requires a previous test output JSON. "
                    "Set --failed_json or ensure --test_output_json exists from a prior run."
                )
            previous = StructuredData(json_path=failed_json)
            failed_ids = previous.failed_ids
            if not failed_ids:
                self._ctx.log("No failed tests in previous run — nothing to re-run.")
                empty = StructuredData(data={"version": "0.1", "tests": []})
                empty.calculate_summary_data()
                return empty
            test_cases = [
                tc for tc in test_cases if hasattr(tc, "_test_id") and f"{tc._test_id}_{tc.index}" in failed_ids
            ]

        test_results = self._collect_test_results(test_cases, test_timeout)
        tests = []
        for test_case, run_response in test_results:
            test_case_data = test_case.structured_test_data(run_response)
            tests.append(test_case_data)
        test_data = {
            "version": "0.1",
            "tests": tests,
        }
        structured_results = StructuredData(data=test_data)
        structured_results.calculate_summary_data()
        return structured_results

    def _collect_test_results(self, test_cases, test_timeout):
        run_responses = self._run_test_cases(test_cases, test_timeout)
        return [(test_case, run_response) for test_case, run_response in zip(test_cases, run_responses)]

    def _run_test_cases(self, test_cases, test_timeout):
        runnables = [test_case.runnable for test_case in test_cases]
        output_collectors = [
            lambda run_response, test_case=test_case: test_case.structured_test_data(run_response)
            for test_case in test_cases
        ]
        with materialized_job_paths(test_cases) as job_paths:
            return self._run(runnables, job_paths, output_collectors, test_timeout=test_timeout)


__all__ = (
    "Engine",
    "BaseEngine",
    "materialized_job_paths",
)
