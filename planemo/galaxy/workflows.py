"""Utilities for Galaxy workflows."""
import json
import os

from collections import namedtuple

import yaml

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


def load_shed_repos(path):
    # TODO: This is crap - doesn't have nested repositories at all.
    if path.endswith(".ga"):
        with open(path, "r") as f:
            workflow = json.load(f)
    else:
        with open(path, "r") as f:
            workflow = yaml.load(f)

    return workflow.get("tools", [])


def install_shed_repos(path, admin_gi):
    tools_info = load_shed_repos(path)
    if tools_info:
        shed_install.install_tools(tools_info, admin_gi, path, default_install_tool_dependencies=False)


def import_workflow(path, admin_gi, user_gi):
    """Import a workflow path to specified Galaxy instance."""
    importer = BioBlendImporterGalaxyInterface(
        admin_gi=admin_gi,
        user_gi=user_gi
    )
    workflow = _raw_dict(path, importer)
    return importer.import_workflow(workflow)


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
