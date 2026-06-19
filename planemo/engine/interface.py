"""Module contianing the :class:`Engine` abstraction."""

import abc
import json
import os
import tempfile
from typing import (
    Callable,
    List,
    Optional,
)

import click

from planemo.exit_codes import EXIT_CODE_UNSUPPORTED_FILE_TYPE
from planemo.io import error
from planemo.runnable import (
    cases,
    RunnableType,
)
from planemo.test.results import StructuredData


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
    def _run(self, runnables, job_path, output_collectors: Optional[List[Callable]] = None):
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
        job_paths = []
        tmp_paths = []
        output_collectors = []
        for test_case in test_cases:
            if test_case.job_path is None:
                job = test_case.job
                with tempfile.NamedTemporaryFile(
                    dir=test_case.tests_directory,
                    suffix=".json",
                    prefix="plnmotmptestjob",
                    delete=False,
                    mode="w+",
                ) as f:
                    tmp_path = f.name
                    job_path = tmp_path
                    tmp_paths.append(tmp_path)
                    json.dump(job, f)
                job_paths.append(job_path)
            else:
                job_paths.append(test_case.job_path)
            output_collectors.append(
                lambda run_response, test_case=test_case: test_case.structured_test_data(run_response)
            )
        try:
            run_responses = self._run(runnables, job_paths, output_collectors)
        finally:
            for tmp_path in tmp_paths:
                os.remove(tmp_path)
        return run_responses


__all__ = (
    "Engine",
    "BaseEngine",
)
