"""Utilities for Galaxy workflows."""
import json
import os

from collections import namedtuple

import yaml

from bioblend.galaxy.client import Client

try:
    from ephemeris import shed_install
except ImportError:
    shed_install = None

try:
    from gxformat2.converter import python_to_workflow
    from gxformat2.interface import BioBlendImporterGalaxyInterface
    from gxformat2.interface import ImporterGalaxyInterface
except ImportError:
    python_to_workflow = None
    BioBlendImporterGalaxyInterface = None
    ImporterGalaxyInterface = object


def load_shed_repos(runnable):
    if runnable.type.name != "galaxy_workflow":
        return []

    # TODO: This is crap - doesn't handle nested repositories at all.
    path = runnable.path
    if path.endswith(".ga"):
        with open(path, "r") as f:
            workflow = json.load(f)
    else:
        with open(path, "r") as f:
            workflow = yaml.load(f)

    return workflow.get("tools", [])


def install_shed_repos(runnable, admin_gi):
    tools_info = load_shed_repos(runnable)
    if tools_info:
        shed_install.install_tools(tools_info, admin_gi, runnable.path, default_install_tool_dependencies=False)


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
