"""Utilities for Galaxy workflows."""
import json
import os
from collections import namedtuple

import yaml
from bioblend.galaxy.client import Client
from ephemeris import generate_tool_list_from_ga_workflow_files
from ephemeris import shed_tools
from gxformat2.converter import python_to_workflow
from gxformat2.interface import BioBlendImporterGalaxyInterface
from gxformat2.interface import ImporterGalaxyInterface


def load_shed_repos(runnable):
    if runnable.type.name != "galaxy_workflow":
        return []

    path = runnable.path
    if path.endswith(".ga"):
        generate_tool_list_from_ga_workflow_files.generate_tool_list_from_workflow([path], "Tools from workflows", "tools.yml")
        with open("tools.yml", "r") as f:
            tools = yaml.load(f)["tools"]

    else:
        # It'd be better to just infer this from the tool shed ID somehow than
        # require explicit annotation like this... I think?
        with open(path, "r") as f:
            workflow = yaml.load(f)

        tools = workflow.get("tools", [])

    return tools


def install_shed_repos(runnable, admin_gi):
    tools_info = load_shed_repos(runnable)
    if tools_info:
        shed_tools._ensure_log_configured("ephemeris")
        install_tool_manager = shed_tools.InstallToolManager(tools_info, admin_gi)
        install_tool_manager.install_repositories()
        if install_tool_manager.errored_repositories:
            raise Exception("Failed to install one or more repositories.")


def import_workflow(path, admin_gi, user_gi, from_path=False):
    """Import a workflow path to specified Galaxy instance."""
    if not from_path:
        importer = BioBlendImporterGalaxyInterface(
            admin_gi=admin_gi,
            user_gi=user_gi
        )
        workflow = _raw_dict(path, importer)
        return importer.import_workflow(workflow)
    else:
        # TODO: Update bioblend to allow from_path.
        payload = dict(
            from_path=path
        )
        workflows_url = user_gi._make_url(user_gi.workflows)
        workflow = Client._post(user_gi.workflows, payload, url=workflows_url)
        return workflow


def _raw_dict(path, importer):
    if path.endswith(".ga"):
        with open(path, "r") as f:
            workflow = json.load(f)
    else:
        workflow_directory = os.path.dirname(path)
        workflow_directory = os.path.abspath(workflow_directory)
        with open(path, "r") as f:
            workflow = yaml.load(f)
            workflow = python_to_workflow(workflow, importer, workflow_directory)

    return workflow


WorkflowOutput = namedtuple("WorkflowOutput", ["order_index", "output_name", "label"])


def describe_outputs(path):
    """Return a list of :class:`WorkflowOutput` objects for target workflow."""
    importer = DummyImporterGalaxyInterface()
    workflow = _raw_dict(path, importer)
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


__all__ = (
    "import_workflow",
    "describe_outputs",
)
