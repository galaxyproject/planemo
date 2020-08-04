"""Describe artifacts that can be run, tested, and linted."""

from __future__ import absolute_import

import abc
import collections
import os
from distutils.dir_util import copy_tree
from enum import auto, Enum

import yaml
from galaxy.tool_util.cwl.parser import workflow_proxy
from galaxy.tool_util.loader_directory import (
    is_a_yaml_with_class,
    looks_like_a_cwl_artifact,
    looks_like_a_data_manager_xml,
    looks_like_a_tool_cwl,
    looks_like_a_tool_xml,
)
from galaxy.tool_util.parser import get_tool_source
from six import (
    add_metaclass,
    python_2_unicode_compatible,
)

from planemo.exit_codes import EXIT_CODE_UNKNOWN_FILE_TYPE, ExitCodeException
from planemo.galaxy.workflows import describe_outputs
from planemo.io import error
from planemo.test import check_output, for_collections

TEST_SUFFIXES = [
    "-tests", "_tests", "-test", "_test"
]
TEST_EXTENSIONS = [".yml", ".yaml", ".json"]

TEST_FILE_NOT_LIST_MESSAGE = ("Invalid test definition file [%s] - file must "
                              "contain a list of tests")
TEST_FIELD_MISSING_MESSAGE = ("Invalid test definition [test #%d in %s] -"
                              "defintion must field [%s].")


class RunnableType(Enum):
    galaxy_tool = auto()
    galaxy_datamanager = auto()
    galaxy_workflow = auto()
    cwl_tool = auto()
    cwl_workflow = auto()
    directory = auto()

    @property
    def has_tools(runnable_type):
        return runnable_type.name in ["galaxy_tool", "galaxy_datamanager", "cwl_tool", "directory"]

    @property
    def is_single_artifact(runnable_type):
        return runnable_type.name not in ["directory"]

    @property
    def test_data_in_parent_dir(runnable_type):
        return runnable_type.name in ["galaxy_datamanager"]

    @property
    def is_galaxy_artifact(runnable_type):
        return "galaxy" in runnable_type.name


_Runnable = collections.namedtuple("Runnable", ["path", "type"])


class Runnable(_Runnable):
    """Abstraction describing tools and workflows."""

    @property
    def test_data_search_path(self):
        """During testing, path to search for test data files."""
        if self.type.name in ['galaxy_datamanager']:
            return os.path.join(os.path.dirname(self.path), os.path.pardir)
        else:
            return self.path

    @property
    def tool_data_search_path(self):
        """During testing, path to search for Galaxy tool data tables."""
        return self.test_data_search_path

    @property
    def data_manager_conf_path(self):
        """Path of a Galaxy data manager configuration for runnable or None."""
        if self.type.name in ['galaxy_datamanager']:
            return os.path.join(os.path.dirname(self.path), os.pardir, 'data_manager_conf.xml')

    @property
    def has_tools(self):
        """Boolean indicating if this runnable corresponds to one or more tools."""
        return _runnable_delegate_attribute('has_tools')

    @property
    def is_single_artifact(self):
        """Boolean indicating if this runnable is a single artifact.

        Currently only directories are considered not a single artifact.
        """
        return _runnable_delegate_attribute('is_single_artifact')


def _runnable_delegate_attribute(attribute):

    @property
    def getter(runnable):
        return getattr(runnable.type, attribute)

    return getter


def _copy_runnable_tree(path, runnable_type, temp_path):
    dir_to_copy = None
    if runnable_type in {RunnableType.galaxy_tool, RunnableType.cwl_tool}:
        dir_to_copy = os.path.dirname(path)
        path = os.path.join(temp_path, os.path.basename(path))
    elif runnable_type == RunnableType.directory:
        dir_to_copy = path
        path = temp_path
    elif runnable_type == RunnableType.galaxy_datamanager:
        dir_to_copy = os.path.join(os.path.dirname(path), os.pardir)
        path_to_data_manager_tool = os.path.relpath(path, dir_to_copy)
        path = os.path.join(temp_path, path_to_data_manager_tool)
    if dir_to_copy:
        copy_tree(dir_to_copy, temp_path, update=True)
    return path


def for_path(path, temp_path=None):
    """Produce a class:`Runnable` for supplied path."""
    runnable_type = None
    if os.path.isdir(path):
        runnable_type = RunnableType.directory
    elif looks_like_a_tool_cwl(path):
        runnable_type = RunnableType.cwl_tool
    elif looks_like_a_data_manager_xml(path):
        runnable_type = RunnableType.galaxy_datamanager
    elif looks_like_a_tool_xml(path):
        runnable_type = RunnableType.galaxy_tool
    elif is_a_yaml_with_class(path, ["GalaxyWorkflow"]):
        runnable_type = RunnableType.galaxy_workflow
    elif path.endswith(".ga"):
        runnable_type = RunnableType.galaxy_workflow
    elif looks_like_a_cwl_artifact(path, ["Workflow"]):
        runnable_type = RunnableType.cwl_workflow
    else:
        # Check to see if it is a Galaxy workflow with a different extension
        try:
            with open(path, "r") as f:
                as_dict = yaml.safe_load(f)
            if as_dict.get("a_galaxy_workflow", False):
                runnable_type = RunnableType.galaxy_workflow
        except Exception:
            pass

    if runnable_type is None:
        error("Unable to determine runnable type for path [%s]" % path)
        raise ExitCodeException(EXIT_CODE_UNKNOWN_FILE_TYPE)

    if temp_path:
        path = _copy_runnable_tree(path, runnable_type, temp_path)

    return Runnable(path, runnable_type)


def for_paths(paths, temp_path=None):
    """Return a specialized list of Runnable objects for paths."""
    return [for_path(path, temp_path=temp_path) for path in paths]


def cases(runnable):
    """Build a `list` of :class:`TestCase` objects for specified runnable."""
    cases = []

    tests_path = _tests_path(runnable)
    if tests_path is None:
        if runnable.type == RunnableType.galaxy_tool:
            tool_source = get_tool_source(runnable.path)
            test_dicts = tool_source.parse_tests_to_dict()
            tool_id = tool_source.parse_id()
            tool_version = tool_source.parse_version()
            for i, test_dict in enumerate(test_dicts.get("tests", [])):
                cases.append(ExternalGalaxyToolTestCase(runnable, tool_id, tool_version, i, test_dict))
        return cases

    tests_directory = os.path.abspath(os.path.dirname(tests_path))

    def normalize_to_tests_path(path):
        if not os.path.isabs(path):
            absolute_path = os.path.join(tests_directory, path)
        else:
            absolute_path = path
        return os.path.normpath(absolute_path)

    with open(tests_path, "r") as f:
        tests_def = yaml.safe_load(f)

    if not isinstance(tests_def, list):
        message = TEST_FILE_NOT_LIST_MESSAGE % tests_path
        raise Exception(message)

    for i, test_def in enumerate(tests_def):
        if "job" not in test_def:
            message = TEST_FIELD_MISSING_MESSAGE % (
                i + 1, tests_path, "job"
            )
            raise Exception(message)
        job_def = test_def["job"]
        if isinstance(job_def, dict):
            job_path = None
            job = job_def
        else:
            job_path = normalize_to_tests_path(job_def)
            job = None

        doc = test_def.get("doc", None)
        output_expectations = test_def.get("outputs", {})
        case = TestCase(
            runnable=runnable,
            tests_directory=tests_directory,
            output_expectations=output_expectations,
            index=i,
            job_path=job_path,
            job=job,
            doc=doc,
        )
        cases.append(case)

    return cases


@add_metaclass(abc.ABCMeta)
class AbstractTestCase(object):
    """Description of a test case for a runnable."""

    def structured_test_data(self, run_response):
        """Result of executing this test case - a "structured_data" dict.

        :rtype: dict
        :return:
                 For example::

                   {
                       "id": "",
                       "has_data": true,
                       "data": {
                           "status": "success", // error, skip,
                           "job": {
                               "command_line": "cat moo",
                               "stdout": "",
                               "stderr": ""
                           },
                           "output_problems": [],
                           "execution_problem": "",
                           "inputs" = {},
                           "problem_log": ""
                       }
                   }
        """


class TestCase(AbstractTestCase):
    """Describe an abstract test case for a specified runnable."""

    def __init__(self, runnable, tests_directory, output_expectations, job_path, job, index, doc):
        """Construct TestCase object from required attributes."""
        self.runnable = runnable
        self.job_path = job_path
        self.job = job
        self.output_expectations = output_expectations
        self.tests_directory = tests_directory
        self.index = index
        self.doc = doc

    def __repr__(self):
        return 'TestCase (%s) for runnable (%s) with job (%s) and expected outputs (%s) in directory (%s) with id (%s)' % \
            (self.doc, self.runnable, self.job, self.output_expectations, self.tests_directory, self.index)

    def structured_test_data(self, run_response):
        """Check a test case against outputs dictionary."""
        output_problems = []
        if run_response.was_successful:
            outputs_dict = run_response.outputs_dict
            execution_problem = None
            for output_id, output_test in self.output_expectations.items():
                if output_id not in outputs_dict:
                    message = "Expected output [%s] not found in results." % output_id
                    output_problems.append(message)
                    continue

                output_value = outputs_dict[output_id]
                output_problems.extend(
                    self._check_output(output_id, output_value, output_test)
                )
            if output_problems:
                status = "failure"
            else:
                status = "success"
        else:
            execution_problem = run_response.error_message
            status = "error"
        data_dict = dict(
            status=status
        )
        if status != "success":
            data_dict["output_problems"] = output_problems
            data_dict["execution_problem"] = execution_problem
        log = run_response.log
        if log is not None:
            data_dict["problem_log"] = log
        job_info = run_response.job_info
        if job_info is not None:
            data_dict["job"] = job_info
        data_dict["inputs"] = self._job
        return dict(
            id=("%s_%s" % (self._test_id, self.index)),
            has_data=True,
            data=data_dict,
        )

    @property
    def _job(self):
        if self.job_path is not None:
            with open(self.job_path, "r") as f:
                return f.read()
        else:
            return self.job

    @property
    def input_ids(self):
        """Labels of inputs specified in test description."""
        return list(self._job.keys())

    @property
    def tested_output_ids(self):
        """Labels of outputs checked in test description."""
        return list(self.output_expectations.keys())

    def _check_output(self, output_id, output_value, output_test):
        output_problems = []
        if not isinstance(output_test, dict):
            if output_test != output_value:
                template = "Output [%s] value [%s] does not match expected value [%s]."
                message = template % (output_id, output_value, output_test)
                output_problems.append(message)
        else:
            if not for_collections(output_test):
                if not isinstance(output_value, dict):
                    message = "Expected file properties for output [%s]" % output_id
                    print(message)
                    print(output_value)
                    output_problems.append(message)
                    return output_problems
                if "path" not in output_value and "location" in output_value:
                    assert output_value["location"].startswith("file://")
                    output_value["path"] = output_value["location"][len("file://"):]
                if "path" not in output_value:
                    message = "No path specified for expected output file [%s]" % output_id
                    output_problems.append(message)
                    print(message)
                    return output_problems
            else:
                output_test["name"] = output_id

            output_problems.extend(
                check_output(
                    self.runnable,
                    output_value,
                    output_test,
                    # TODO: needs kwds in here...
                )
            )

        return output_problems

    @property
    def _test_id(self):
        if self.runnable.type in [
            RunnableType.cwl_tool,
            RunnableType.galaxy_tool,
        ]:
            return get_tool_source(self.runnable.path).parse_id()
        else:
            return os.path.basename(self.runnable.path)


class ExternalGalaxyToolTestCase(AbstractTestCase):
    """Special class of AbstractCase that doesn't use job_path but uses test data from a Galaxy server."""

    def __init__(self, runnable, tool_id, tool_version, test_index, test_dict):
        """Construct TestCase object from required attributes."""
        self.runnable = runnable
        self.tool_id = tool_id
        self.tool_version = tool_version
        self.test_index = test_index
        self.test_dict = test_dict

    def structured_test_data(self, run_response):
        """Just return the structured_test_data generated from galaxy-tool-util for this test variant."""
        return run_response


def _tests_path(runnable):
    if not runnable.is_single_artifact:
        raise NotImplementedError("Tests for directories are not yet implemented.")

    runnable_path = runnable.path
    base, _ = os.path.splitext(runnable_path)

    for test_suffix in TEST_SUFFIXES:
        for test_extension in TEST_EXTENSIONS:
            test_path = base + test_suffix + test_extension
            if os.path.exists(test_path):
                return test_path

    return None


def get_outputs(runnable):
    """Return a list of :class:`RunnableOutput` objects for this runnable."""
    if not runnable.is_single_artifact:
        raise NotImplementedError("Cannot generate outputs for a directory.")
    if runnable.type in [RunnableType.galaxy_tool, RunnableType.cwl_tool]:
        tool_source = get_tool_source(runnable.path)
        # TODO: do something with collections at some point
        output_datasets, _ = tool_source.parse_outputs(None)
        outputs = [ToolOutput(o) for o in output_datasets.values()]
        return outputs
    elif runnable.type == RunnableType.galaxy_workflow:
        workflow_outputs = describe_outputs(runnable.path)
        return [GalaxyWorkflowOutput(o) for o in workflow_outputs]
    elif runnable.type == RunnableType.cwl_workflow:
        workflow = workflow_proxy(runnable.path, strict_cwl_validation=False)
        return [CwlWorkflowOutput(label) for label in workflow.output_labels]
    else:
        raise NotImplementedError("Getting outputs for this artifact type is not yet supported.")


@add_metaclass(abc.ABCMeta)
class RunnableOutput(object):
    """Description of a single output of an execution of a Runnable."""

    @abc.abstractproperty
    def get_id(self):
        """An identifier that describes this output."""


class ToolOutput(RunnableOutput):
    """Implementation of RunnableOutput corresponding to Galaxy tool outputs."""

    def __init__(self, tool_output):
        self._tool_output = tool_output

    def get_id(self):
        return self._tool_output.name


class GalaxyWorkflowOutput(RunnableOutput):
    """Implementation of RunnableOutput corresponding to Galaxy workflow outputs."""

    def __init__(self, workflow_output):
        self._workflow_output = workflow_output

    def get_id(self):
        return self._workflow_output.label

    @property
    def workflow_output(self):
        return self._workflow_output


class CwlWorkflowOutput(RunnableOutput):
    """Implementation of RunnableOutput corresponding to CWL outputs."""

    def __init__(self, label):
        self._label = label

    def get_id(self):
        return self._label


@add_metaclass(abc.ABCMeta)
class RunResponse(object):
    """Description of an attempt for an engine to execute a Runnable."""

    @abc.abstractproperty
    def was_successful(self):
        """Indicate whether an error was encountered while executing this runnable.

        If successful, response should conform to the SuccessfulRunResponse interface,
        otherwise it will conform to the ErrorRunResponse interface.
        """

    @abc.abstractproperty
    def job_info(self):
        """If job information is available, return as dictionary."""

    @abc.abstractproperty
    def log(self):
        """If engine related log is available, return as text data."""


@add_metaclass(abc.ABCMeta)
class SuccessfulRunResponse(RunResponse):
    """Description of the results of an engine executing a Runnable."""

    def was_successful(self):
        """Return `True` to indicate this run was successful."""
        return True

    @abc.abstractproperty
    def outputs_dict(self):
        """Return a dict of output descriptions."""


@python_2_unicode_compatible
class ErrorRunResponse(RunResponse):
    """Description of an error while attempting to execute a Runnable."""

    def __init__(self, error_message, job_info=None, log=None):
        """Create an ErrorRunResponse with specified error message."""
        self._error_message = error_message
        self._job_info = job_info
        self._log = log

    @property
    def error_message(self):
        """Error message describing the problem with execution of the runnable."""
        return self._error_message

    @property
    def was_successful(self):
        """Return `False` to indicate this run was successful."""
        return False

    @property
    def job_info(self):
        """Return potentially null stored `job_info` dict."""
        return self._job_info

    @property
    def log(self):
        """Return potentially null stored `log` text."""
        return self._log

    def __str__(self):
        """Print a helpful error description of run."""
        message = "Run failed with message [%s]" % self.error_message
        log = self.log
        if log:
            message += " and log [%s]" % log
        return message


__all__ = (
    "cases",
    "ErrorRunResponse",
    "for_path",
    "for_paths",
    "get_outputs",
    "Runnable",
    "RunnableType",
    "RunResponse",
    "RunnableOutput",
    "SuccessfulRunResponse",
    "TestCase",
)
