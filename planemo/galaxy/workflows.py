"""Utilities for Galaxy workflows."""
import json
import os
from collections import namedtuple
from urllib.parse import urlparse

import yaml
from ephemeris import generate_tool_list_from_ga_workflow_files
from ephemeris import shed_tools
from gxformat2.converter import python_to_workflow
from gxformat2.interface import BioBlendImporterGalaxyInterface
from gxformat2.interface import ImporterGalaxyInterface
from gxformat2.normalize import inputs_normalized, outputs_normalized

from planemo.io import warn

FAILED_REPOSITORIES_MESSAGE = "Failed to install one or more repositories."
GALAXY_WORKFLOWS_PREFIX = "gxid://workflows/"


def load_shed_repos(runnable):
    if runnable.type.name != "galaxy_workflow":
        return []
    path = runnable.path
    if path.endswith(".ga"):
        generate_tool_list_from_ga_workflow_files.generate_tool_list_from_workflow([path], "Tools from workflows", "tools.yml")
        with open("tools.yml", "r") as f:
            tools = yaml.safe_load(f)["tools"]

    else:
        # It'd be better to just infer this from the tool shed ID somehow than
        # require explicit annotation like this... I think?
        with open(path, "r") as f:
            workflow = yaml.safe_load(f)

        tools = workflow.get("tools", [])

    return tools


def install_shed_repos(runnable, admin_gi,
                       ignore_dependency_problems,
                       install_tool_dependencies=False,
                       install_resolver_dependencies=True,
                       install_repository_dependencies=True):
    tools_info = load_shed_repos(runnable)
    if tools_info:
        install_tool_manager = shed_tools.InstallRepositoryManager(admin_gi)
        install_results = install_tool_manager.install_repositories(tools_info,
                                                                    default_install_tool_dependencies=install_tool_dependencies,
                                                                    default_install_resolver_dependencies=install_resolver_dependencies,
                                                                    default_install_repository_dependencies=install_repository_dependencies)
        if install_results.errored_repositories:
            if ignore_dependency_problems:
                warn(FAILED_REPOSITORIES_MESSAGE)
            else:
                raise Exception(FAILED_REPOSITORIES_MESSAGE)


def import_workflow(path, admin_gi, user_gi, from_path=False):
    """Import a workflow path to specified Galaxy instance."""
    if not from_path:
        importer = BioBlendImporterGalaxyInterface(
            admin_gi=admin_gi,
            user_gi=user_gi
        )
        workflow = _raw_dict(path, importer)
        return user_gi.workflows.import_workflow_dict(workflow)
    else:
        path = os.path.abspath(path)
        workflow = user_gi.workflows.import_workflow_from_local_path(path)
        return workflow


def _raw_dict(path, importer=None):
    if path.endswith(".ga"):
        with open(path, "r") as f:
            workflow = json.load(f)
    else:
        if importer is None:
            importer = DummyImporterGalaxyInterface()

        workflow_directory = os.path.dirname(path)
        workflow_directory = os.path.abspath(workflow_directory)
        with open(path, "r") as f:
            workflow = yaml.safe_load(f)
            workflow = python_to_workflow(workflow, importer, workflow_directory)

    return workflow


def find_tool_ids(path):
    tool_ids = set()
    workflow = _raw_dict(path)

    def register_tool_ids(tool_ids, workflow):
        for step in workflow["steps"].values():
            if step.get('subworkflow'):
                register_tool_ids(tool_ids, step['subworkflow'])
            elif step.get("tool_id"):
                tool_ids.add(step['tool_id'])

    register_tool_ids(tool_ids, workflow)

    return list(tool_ids)


WorkflowOutput = namedtuple("WorkflowOutput", ["order_index", "output_name", "label"])


def remote_runnable_to_workflow_id(runnable):
    assert runnable.is_remote_workflow_uri
    parse_result = urlparse(runnable.uri)
    return parse_result.path[1:]


def describe_outputs(runnable, gi=None):
    """Return a list of :class:`WorkflowOutput` objects for target workflow."""
    if runnable.uri.startswith(GALAXY_WORKFLOWS_PREFIX):
        workflow_id = remote_runnable_to_workflow_id(runnable)
        assert gi is not None
        workflow = get_dict_from_workflow(gi, workflow_id)
    else:
        workflow = _raw_dict(runnable.path)

    outputs = []
    for (order_index, step) in workflow["steps"].items():
        step_outputs = step.get("workflow_outputs", [])
        for step_output in step_outputs:
            output = WorkflowOutput(
                int(order_index),
                step_output["output_name"],
                step_output["label"],
            )
            outputs.append(output)
    return outputs


class DummyImporterGalaxyInterface(ImporterGalaxyInterface):

    def import_workflow(self, workflow, **kwds):
        return None


def input_labels(workflow_path):
    """Get normalized labels for workflow artifact regardless of format."""
    steps = inputs_normalized(workflow_path=workflow_path)
    labels = []
    for step in steps:
        step_id = input_label(step)
        if step_id:
            labels.append(step_id)
    return labels


def required_input_steps(workflow_path):
    try:
        steps = inputs_normalized(workflow_path=workflow_path)
    except Exception:
        raise Exception("Input workflow could not be successfully normalized - try linting with planemo workflow_lint.")
    required_steps = []
    for input_step in steps:
        if input_step.get("optional", False) or input_step.get("default"):
            continue
        required_steps.append(input_step)
    return required_steps


def required_input_labels(workflow_path):
    return map(input_label, required_input_steps(workflow_path))


def input_label(input_step):
    """Get the normalized label of a step returned from inputs_normalized."""
    step_id = input_step.get("id") or input_step.get("label")
    return step_id


def output_labels(workflow_path):
    outputs = outputs_normalized(workflow_path=workflow_path)
    return [o["id"] for o in outputs]


def output_stubs_for_workflow(workflow_path):
    """
    Return output labels and class.
    """
    outputs = {}
    for label in output_labels(workflow_path):
        if not label.startswith('_anonymous_'):
            outputs[label] = {'class': ''}
    return outputs


def job_template(workflow_path):
    """Return a job template for specified workflow.

    A dictionary describing non-optional inputs that must be specified to
    run the workflow.
    """
    template = {}
    for required_input_step in required_input_steps(workflow_path):
        i_label = input_label(required_input_step)
        input_type = required_input_step["type"]
        if input_type == "data":
            template[i_label] = {
                "class": "File",
                "path": "todo_test_data_path.ext",
            }
        elif input_type == "collection":
            template[i_label] = {
                "class": "Collection",
                "collection_type": "list",
                "elements": [
                    {
                        "class": "File",
                        "identifier": "todo_element_name",
                        "path": "todo_test_data_path.ext",
                    }
                ],
            }
        elif input_type in ['string', 'int', 'float', 'boolean', 'color']:
            template[i_label] = "todo_param_value"
        else:
            template[i_label] = {
                "TODO",  # Does this work yet?
            }
    return template


def new_workflow_associated_path(workflow_path, suffix="tests"):
    """Generate path for test or job YAML file next to workflow."""
    base, input_ext = os.path.splitext(workflow_path)
    # prefer -tests.yml but if the author uses underscores or .yaml respect that.
    sep = "-"
    if "_" in base and "-" not in base:
        sep = "_"
    ext = "yml"
    if "yaml" in input_ext:
        ext = "yaml"
    return base + sep + suffix + "." + ext


def get_dict_from_workflow(gi, workflow_id):
    return gi.workflows.export_workflow_dict(workflow_id)


def rewrite_job_file(input_file, output_file, job):
    """Rewrite a job file with galaxy_ids for upload_data subcommand"""
    with open(input_file, "r") as f:
        job_contents = yaml.safe_load(f)
        for job_input, job_input_name in job_contents.items():
            if type(job[job_input]) == dict:  # dataset or collection
                job_contents[job_input] = {
                    "class": job_input_name["class"],
                    "galaxy_id": job[job_input]["id"]
                }
            # else: presumably a parameter, no need to modify
    with open(output_file, "w") as f:
        yaml.dump(job_contents, f)


__all__ = (
    "import_workflow",
    "describe_outputs",
)
