""" Planemo specific utilities for dealing with conda, extending Galaxy's
features with planemo specific idioms.
"""
from __future__ import absolute_import

from galaxy.tools.deps import conda_util
from planemo.io import shell

from galaxy.tools.deps.requirements import parse_requirements_from_xml
from galaxy.tools.loader_directory import load_tool_elements_from_path


def build_conda_context(**kwds):
    """ Build a Galaxy CondaContext tailored to planemo use
    and common command-line arguments.
    """
    conda_prefix = kwds.get("conda_prefix", None)
    use_planemo_shell = kwds.get("use_planemo_shell_exec", True)
    ensure_channels = kwds.get("conda_ensure_channels", "")
    shell_exec = shell if use_planemo_shell else None
    return conda_util.CondaContext(conda_prefix=conda_prefix,
                                   ensure_channels=ensure_channels,
                                   shell_exec=shell_exec)


def collect_conda_targets(path, found_tool_callback=None, conda_context=None):
    conda_targets = []
    for (tool_path, tool_xml) in load_tool_elements_from_path(path):
        if found_tool_callback:
            found_tool_callback(tool_path)
        requirements, containers = parse_requirements_from_xml(tool_xml)
        conda_targets.extend(conda_util.requirements_to_conda_targets(requirements))
    return conda_targets
