"""Module contianing the :class:`Engine` abstraction."""

import abc
import json
import os
import tempfile
from typing import List

from six import add_metaclass

from planemo.exit_codes import EXIT_CODE_UNSUPPORTED_FILE_TYPE
from planemo.io import error
from planemo.runnable import (
    cases,
    for_path,
    RunnableType,
)
from planemo.test.results import StructuredData


@add_metaclass(abc.ABCMeta)
class Engine(object):
    """Abstract description of an external process for running tools or workflows.
    """

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

    handled_runnable_types = []  # type: List[RunnableType]

    def __init__(self, ctx, **kwds):
        """Store context and kwds."""
        self._ctx = ctx
        self._kwds = kwds

    def can_run(self, runnable):
        """Use subclass's ``handled_runnable_types`` variable to infer ``can_run``."""
        return runnable.type in self.handled_runnable_types

    def cleanup(self):
        """Default no-op cleanup method."""

    def run(self, path, job_path):
        """Run a job using a compatible artifact (workflow or tool)."""
        runnable = for_path(path)
        self._check_can_run(runnable)
        run_response = self._run(runnable, job_path)
        return run_response

    @abc.abstractmethod
    def _run(self, runnable, job_path):
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

    def test(self, runnables):
        """Test runnable artifacts (workflow or tool)."""
        self._check_can_run_all(runnables)
        test_cases = [t for tl in map(cases, runnables) for t in tl]
        test_results = self._collect_test_results(test_cases)
        tests = []
        for (test_case, run_response) in test_results:
            test_case_data = test_case.structured_test_data(run_response)
            tests.append(test_case_data)
        test_data = {
            'version': '0.1',
            'tests': tests,
        }
        structured_results = StructuredData(data=test_data)
        structured_results.calculate_summary_data()
        return structured_results

    def _collect_test_results(self, test_cases):
        test_results = []
        for test_case in test_cases:
            self._ctx.vlog(
                "Running tests %s" % test_case
            )
            run_response = self._run_test_case(test_case)
            self._ctx.vlog(
                "Test case [%s] resulted in run response [%s]",
                test_case,
                run_response,
            )
            test_results.append((test_case, run_response))
        return test_results

    def _run_test_case(self, test_case):
        runnable = test_case.runnable
        job_path = test_case.job_path
        tmp_path = None
        if job_path is None:
            job = test_case.job
            f = tempfile.NamedTemporaryFile(
                dir=test_case.tests_directory,
                suffix=".json",
                prefix="plnmotmptestjob",
                delete=False,
                mode="w+",
            )
            tmp_path = f.name
            job_path = tmp_path
            json.dump(job, f)
            f.close()
        try:
            run_response = self._run(runnable, job_path)
        finally:
            if tmp_path:
                os.remove(tmp_path)
        return run_response


__all__ = (
    "Engine",
    "BaseEngine",
)
