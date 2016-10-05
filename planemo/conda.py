"""Planemo specific utilities for dealing with conda.

The extend Galaxy/galaxy-lib's features with planemo specific idioms.
"""
from __future__ import absolute_import

import os

from galaxy.tools.deps import conda_util

from planemo.io import shell
from planemo.tools import yield_tool_sources_on_paths


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


def collect_conda_targets(ctx, paths, found_tool_callback=None, conda_context=None):
    """Load CondaTarget objects from supplied artifact sources.

    If a tool contains more than one requirement, the requirements will each
    appear once in the output.
    """
    conda_targets = set([])
    for (tool_path, tool_source) in yield_tool_sources_on_paths(ctx, paths):
        if found_tool_callback:
            found_tool_callback(tool_path)
        for target in tool_source_conda_targets(tool_source):
            conda_targets.add(target)
    return conda_targets


def collect_conda_target_lists(ctx, paths, found_tool_callback=None):
    """Load CondaTarget lists from supplied artifact sources.

    If a tool contains more than one requirement, the requirements will all
    appear together as one list element of the output list.
    """
    conda_target_lists = set([])
    for (tool_path, tool_source) in yield_tool_sources_on_paths(ctx, paths):
        if found_tool_callback:
            found_tool_callback(tool_path)
        conda_target_lists.add(frozenset(tool_source_conda_targets(tool_source)))
    return conda_target_lists


def tool_source_conda_targets(tool_source):
    """Load CondaTarget object from supplied abstract tool source."""
    requirements, _ = tool_source.parse_requirements_and_containers()
    return conda_util.requirements_to_conda_targets(requirements)


__all__ = [
    "build_conda_context",
    "collect_conda_targets",
    "collect_conda_target_lists",
    "tool_source_conda_targets",
]
