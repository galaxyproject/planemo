"""Planemo specific utilities for dealing with conda.

The extend Galaxy/galaxy-lib's features with planemo specific idioms.
"""
from __future__ import absolute_import

import os

from galaxy.tools.deps import conda_util
from galaxy.tools.deps.requirements import parse_requirements_from_xml
from galaxy.tools.loader_directory import load_tool_elements_from_path

from planemo.io import shell


def build_conda_context(ctx, **kwds):
    """Build a galaxy-lib CondaContext tailored to planemo use.

    Using planemo's common command-line/globa config options.
    """
    condarc_override_default = os.path.join(ctx.workspace, "condarc")
    conda_prefix = kwds.get("conda_prefix", None)
    use_planemo_shell = kwds.get("use_planemo_shell_exec", True)
    ensure_channels = kwds.get("conda_ensure_channels", "")
    condarc_override = kwds.get("condarc", condarc_override_default)
    shell_exec = shell if use_planemo_shell else None
    return conda_util.CondaContext(conda_prefix=conda_prefix,
                                   ensure_channels=ensure_channels,
                                   condarc_override=condarc_override,
                                   shell_exec=shell_exec)


def collect_conda_targets(path, found_tool_callback=None, conda_context=None):
    """Load CondaTarget objects from supplied artifact sources."""
    conda_targets = []
    for (tool_path, tool_xml) in load_tool_elements_from_path(path):
        if found_tool_callback:
            found_tool_callback(tool_path)
        requirements, containers = parse_requirements_from_xml(tool_xml)
        conda_targets.extend(conda_util.requirements_to_conda_targets(requirements))
    return conda_targets
