import inspect
import json
import os
import re
from collections import OrderedDict
from typing import (
    Any,
    Dict,
    Iterable,
    Iterator,
    List,
    Optional,
    TYPE_CHECKING,
    Union,
)

import yaml
from galaxy.tool_util.lint import LintContext
from galaxy.tool_util.loader_directory import EXCLUDE_WALK_DIRS
from galaxy.tool_util.verify import asserts
from gxformat2.lint import (
    lint_format2,
    lint_ga,
)
from gxformat2.yaml import ordered_load

from planemo.exit_codes import (
    EXIT_CODE_GENERIC_FAILURE,
    EXIT_CODE_OK,
)
from planemo.galaxy.workflows import (
    input_labels,
    output_labels,
    required_input_labels,
)
from planemo.runnable import (
    cases,
    for_path,
    TestCase,
)
from planemo.shed import DOCKSTORE_REGISTRY_CONF

if TYPE_CHECKING:
    from planemo.cli import PlanemoCliContext

POTENTIAL_WORKFLOW_FILES = re.compile(r"^.*(\.yml|\.yaml|\.ga)$")
DOCKSTORE_REGISTRY_CONF_VERSION = "1.2"


class WorkflowLintContext(LintContext):
    # Setup training topic for linting - probably should pass this through
    # from click arguments.
    training_topic = None


def generate_dockstore_yaml(directory: str, publish: bool = True) -> str:
    workflows = []
    for workflow_path in find_workflow_descriptions(directory):
        test_parameter_path = f"{workflow_path.rsplit('.', 1)[0]}-tests.yml"
        workflow_entry: Dict[str, Any] = {
            # TODO: support CWL
            "subclass": "Galaxy",
            "publish": publish,
            "name": "main",
            "primaryDescriptorPath": f"/{os.path.relpath(workflow_path, directory)}",
        }
        if os.path.exists(test_parameter_path):
            workflow_entry["testParameterFiles"] = [f"/{os.path.relpath(test_parameter_path, directory)}"]
        workflows.append(workflow_entry)
    # Force version to the top of file but serializing rest of config separately
    contents = "version: %s\n" % DOCKSTORE_REGISTRY_CONF_VERSION
    contents += yaml.dump({"workflows": workflows})
    return contents


def lint_workflow_artifacts_on_paths(
    ctx: "PlanemoCliContext", paths: Iterable[str], lint_args: Dict[str, Union[str, List[str]]]
) -> int:
    report_level = lint_args["level"]
    lint_context = WorkflowLintContext(report_level, skip_types=lint_args["skip_types"])
    for path in paths:
        _lint_workflow_artifacts_on_path(lint_context, path, lint_args)

    if lint_context.failed(lint_args["fail_level"]):
        return EXIT_CODE_GENERIC_FAILURE
    else:
        return EXIT_CODE_OK


def _lint_workflow_artifacts_on_path(
    lint_context: WorkflowLintContext, path: str, lint_args: Dict[str, Union[str, List[str]]]
) -> None:
    for potential_workflow_artifact_path in find_potential_workflow_files(path):
        if os.path.basename(potential_workflow_artifact_path) == DOCKSTORE_REGISTRY_CONF:
            lint_context.lint("lint_dockstore", _lint_dockstore_config, potential_workflow_artifact_path)

        elif looks_like_a_workflow(potential_workflow_artifact_path):

            def structure(path, lint_context):
                with open(path) as f:
                    workflow_dict = ordered_load(f)
                workflow_class = workflow_dict.get("class")
                lint_func = lint_format2 if workflow_class == "GalaxyWorkflow" else lint_ga
                lint_func(lint_context, workflow_dict, path=path)

            lint_context.lint("lint_structure", structure, potential_workflow_artifact_path)
            lint_context.lint("lint_best_practices", _lint_best_practices, potential_workflow_artifact_path)
            lint_context.lint("lint_tests", _lint_tsts, potential_workflow_artifact_path)
        else:
            # Allow linting ro crates and such also
            pass


# misspell for pytest
def _lint_tsts(path: str, lint_context: WorkflowLintContext) -> None:
    runnables = for_path(path, return_all=True)
    if not isinstance(runnables, list):
        runnables = [runnables]
    for runnable in runnables:
        test_cases = cases(runnable)
        if len(test_cases) == 0:
            lint_context.warn("Workflow missing test cases.")
            return
        all_tests_valid = True
        for test_case in test_cases:
            if isinstance(test_case, TestCase):
                if not _lint_case(path, test_case, lint_context):
                    all_tests_valid = False
            else:
                lint_context.warn(f"Test case of type {type(test_case)} not currently supported.")
                all_tests_valid = False

        if all_tests_valid:
            lint_context.valid(f"Tests appear structurally correct for {runnable.path}")


def _lint_best_practices(path: str, lint_context: WorkflowLintContext) -> None:  # noqa: C901
    """
    This function duplicates the checks made by Galaxy's best practices panel:
    https://github.com/galaxyproject/galaxy/blob/5396bb15fe8cfcf2e89d46c1d061c49b60e2f0b1/client/src/components/Workflow/Editor/Lint.vue
    """

    def check_json_for_untyped_params(j):
        values = j if type(j) == list else j.values()
        for value in values:
            if type(value) in [list, dict, OrderedDict]:
                if check_json_for_untyped_params(value):
                    return True
            elif type(value) == str:
                if re.match(r"\$\{.+?\}", value):
                    return True
        return False

    with open(path) as f:
        workflow_dict = ordered_load(f)

    steps = workflow_dict.get("steps", {})

    # annotation
    if not workflow_dict.get("annotation"):
        lint_context.warn("Workflow is not annotated.")

    # creator
    if not len(workflow_dict.get("creator", [])) > 0:
        lint_context.warn("Workflow does not specify a creator.")

    # license
    if not workflow_dict.get("license"):
        lint_context.warn("Workflow does not specify a license.")

    # checks on individual steps
    for step in steps.values():
        print(step)
        # disconnected inputs
        for input in step.get("inputs", []):
            if input.get("name") not in step.get("input_connections"):  # TODO: check optional
                lint_context.warn(
                    f"Input {input.get('name')} of workflow step {step.get('annotation') or step.get('id')} is disconnected."
                )

        # missing metadata
        if not step.get("annotation"):
            lint_context.warn(f"Workflow step with ID {step.get('id')} has no annotation.")
        if not step.get("label"):
            lint_context.warn(f"Workflow step with ID {step.get('id')} has no label.")

        # untyped parameters
        if workflow_dict.get("class") == "GalaxyWorkflow":
            tool_state = step.get("tool_state", {})
            pjas = step.get("out", {})
        else:
            tool_state = json.loads(step.get("tool_state", "{}"))
            pjas = step.get("post_job_actions", {})

        if check_json_for_untyped_params(tool_state):
            lint_context.warn(f"Workflow step with ID {step.get('id')} specifies an untyped parameter as an input.")

        if check_json_for_untyped_params(pjas):
            lint_context.warn(
                f"Workflow step with ID {step.get('id')} specifies an untyped parameter in the post-job actions."
            )

        # unlabeled outputs are checked by gxformat2, no need to check here


def _lint_case(path: str, test_case: TestCase, lint_context: WorkflowLintContext) -> bool:
    test_valid = True

    i_labels = input_labels(workflow_path=path)
    job_keys = test_case.input_ids
    for key in job_keys:
        if key not in i_labels:
            # consider an error instead?
            lint_context.warn(
                f"Unknown workflow input in test job definition [{key}], workflow inputs are [{i_labels}]"
            )
            test_valid = False

    # check non-optional parameters are set
    for required_label in required_input_labels(path):
        if required_label not in job_keys:
            template = "Non-optional input has no value specified in workflow test job [%s], job specifies inputs [%s]"
            lint_context.error(template % (required_label, job_keys))
            test_valid = False

    for input_id, input_def in test_case._job.items():
        if not _tst_input_valid(test_case, input_id, input_def, lint_context):
            test_valid = False

    test_output_ids = test_case.tested_output_ids
    o_labels = output_labels(path)
    found_valid_expectation = False
    for test_output_id in test_output_ids:
        if test_output_id not in o_labels:
            template = "Test found for unknown workflow output [%s], workflow outputs [%s]"
            lint_context.error(template % (test_output_id, o_labels))
            test_valid = False
        else:
            found_valid_expectation = True
        # TODO: validate structure of test expectations

        assertion_definitions = test_case.output_expectations[test_output_id].get("asserts")
        if not _check_test_assertions(lint_context, assertion_definitions):
            test_valid = False

    if not found_valid_expectation:
        lint_context.warn("Found no valid test expectations for workflow test")
        test_valid = False

    return test_valid


def _check_test_assertions(
    lint_context: WorkflowLintContext, assertion_definitions: Optional[Dict[str, Dict[str, Any]]]
) -> bool:
    # we are already in Python, not XML, so it is simpler to lint assertions by checking against the
    # Python functions directly, rather than checking against galaxy.xsd as for tool linting
    assertions_valid = True
    if assertion_definitions:
        for assertion_name, assertion_params in assertion_definitions.items():
            function = asserts.assertion_functions[f"assert_{assertion_name}"]
            signature = inspect.signature(function)
            try:
                # try mapping the function with the attributes supplied and check for TypeError
                signature.bind("", **assertion_params)
            except AssertionError:
                pass
            except TypeError as e:
                lint_context.error(f"Invalid assertion in tests: assert_{assertion_name} {str(e)}")
                assertions_valid = False
    return assertions_valid


def _tst_input_valid(
    test_case: TestCase,
    input_id: str,
    input_def: Dict[str, Any],
    lint_context: WorkflowLintContext,
) -> bool:
    if type(input_def) == dict:  # else assume it is a parameter
        clazz = input_def.get("class")
        if clazz == "File":
            input_path = input_def.get("path")
            if input_path:
                if not os.path.isabs(input_path):
                    input_path = os.path.join(test_case.tests_directory, input_path)
                if not os.path.exists(input_path):
                    message = "Test referenced File path [%s] not found" % input_path
                    lint_context.warn(message)
                    return False
        elif clazz == "Collection":
            for elem in input_def.get("elements", []):
                elem_valid = _tst_input_valid(test_case, input_id, elem, lint_context)
                if not elem_valid:
                    return False
    return True


def _lint_dockstore_config(path: str, lint_context: WorkflowLintContext) -> None:
    dockstore_yaml = None
    try:
        with open(path) as f:
            dockstore_yaml = yaml.safe_load(f)
    except Exception:
        lint_context.error("Invalid YAML found in %s" % DOCKSTORE_REGISTRY_CONF)
        return

    if not isinstance(dockstore_yaml, dict):
        lint_context.error("Invalid YAML contents found in %s" % DOCKSTORE_REGISTRY_CONF)
        return

    if "workflows" not in dockstore_yaml:
        lint_context.error("Invalid YAML contents found in %s, no workflows defined" % DOCKSTORE_REGISTRY_CONF)
        return

    workflow_entries = dockstore_yaml.get("workflows")
    if not isinstance(workflow_entries, list):
        lint_context.error("Invalid YAML contents found in %s, workflows not a list" % DOCKSTORE_REGISTRY_CONF)
        return

    for workflow_entry in workflow_entries:
        _lint_dockstore_workflow_entry(lint_context, os.path.dirname(path), workflow_entry)


def _lint_dockstore_workflow_entry(
    lint_context: WorkflowLintContext, directory: str, workflow_entry: Dict[str, Any]
) -> None:
    if not isinstance(workflow_entry, dict):
        lint_context.error("Invalid YAML contents found in %s, workflow entry not a dict" % DOCKSTORE_REGISTRY_CONF)
        return

    found_errors = False
    for required_key in ["primaryDescriptorPath", "subclass"]:
        if required_key not in workflow_entry:
            lint_context.error(f"{DOCKSTORE_REGISTRY_CONF} workflow entry missing required key {required_key}")
            found_errors = True

    for recommended_key in ["testParameterFiles"]:
        if recommended_key not in workflow_entry:
            lint_context.warn(f"{DOCKSTORE_REGISTRY_CONF} workflow entry missing recommended key {recommended_key}")

    if found_errors:
        # Don't do the rest of the validation for a broken file.
        return

    # TODO: validate subclass
    descriptor_path = workflow_entry["primaryDescriptorPath"]
    test_files = workflow_entry.get("testParameterFiles", [])

    for referenced_file in [descriptor_path] + test_files:
        referenced_path = os.path.join(directory, referenced_file[1:])
        if not os.path.exists(referenced_path):
            lint_context.error(f"{DOCKSTORE_REGISTRY_CONF} workflow entry references absent file {referenced_file}")


def looks_like_a_workflow(path: str) -> bool:
    """Return boolean indicating if this path looks like a workflow."""
    if POTENTIAL_WORKFLOW_FILES.match(os.path.basename(path)):
        with open(path) as f:
            workflow_dict = ordered_load(f)
        if not isinstance(workflow_dict, dict):
            # Not exactly right - could have a #main def - do better and sync with Galaxy.
            return False
        return bool(workflow_dict.get("class") == "GalaxyWorkflow" or workflow_dict.get("a_galaxy_workflow"))
    return False


def find_workflow_descriptions(directory: str) -> Iterator[str]:
    for potential_workflow_artifact_path in find_potential_workflow_files(directory):
        if looks_like_a_workflow(potential_workflow_artifact_path):
            yield potential_workflow_artifact_path


def find_potential_workflow_files(directory: str) -> List[str]:
    """Return a list of potential workflow files in a directory."""
    if not os.path.exists(directory):
        raise ValueError(f"Directory not found {directory}")

    matches = []
    if os.path.isdir(directory):
        for root, dirnames, filenames in os.walk(directory):
            # exclude some directories (like .hg) from traversing
            dirnames[:] = [dir for dir in dirnames if dir not in EXCLUDE_WALK_DIRS]
            for filename in filenames:
                if POTENTIAL_WORKFLOW_FILES.match(filename):
                    matches.append(os.path.join(root, filename))
    else:
        matches.append(directory)
    return matches
