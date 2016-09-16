"""Planemo-specific wrappers around galaxy-lib tool functionality."""
from __future__ import absolute_import

import os
import sys
import traceback

from galaxy.tools import loader_directory

from planemo.io import error, info

is_tool_load_error = loader_directory.is_tool_load_error
SKIP_XML_MESSAGE = "Skipping XML file - does not appear to be a tool %s."
SHED_FILES = ["tool_dependencies.xml", "repository_dependencies.xml"]


def yield_tool_sources_on_paths(ctx, paths, recursive=False):
    for path in paths:
        for (tool_path, tool_source) in yield_tool_sources(ctx, path, recursive):
            yield (tool_path, tool_source)


def yield_tool_sources(ctx, path, recursive=False):
    tools = load_tool_sources_from_path(
        path,
        recursive,
        register_load_errors=True,
    )
    for (tool_path, tool_source) in tools:
        if is_tool_load_error(tool_source):
            yield (tool_path, tool_source)
            continue
        if not _is_tool_source(ctx, tool_path, tool_source):
            continue
        yield (tool_path, tool_source)


def load_tool_sources_from_path(path, recursive, register_load_errors=False):
    """Generator for tool sources on a path."""
    return loader_directory.load_tool_sources_from_path(
        path,
        _load_exception_handler,
        recursive=recursive,
        register_load_errors=register_load_errors,
    )


# TODO: replace usage of this with tool source so planemo
# defers XML details to galaxy-lib.
def load_tool_elements_from_path(path, recursive, register_load_errors=False):
    """Generator for tool XML elements on a path."""
    return loader_directory.load_tool_elements_from_path(
        path,
        _load_exception_handler,
        recursive=recursive,
        register_load_errors=register_load_errors,
    )


def _load_exception_handler(path, exc_info):
    error("Error loading tool with path %s" % path)
    traceback.print_exception(*exc_info, limit=1, file=sys.stderr)


def _is_tool_source(ctx, tool_path, tool_source):
    if os.path.basename(tool_path) in SHED_FILES:
        return False
    root = getattr(tool_source, "root", None)
    if root is not None:
        if root.tag != "tool":
            if ctx.verbose:
                info(SKIP_XML_MESSAGE % tool_path)
            return False
    return True


__all__ = [
    "load_tool_elements_from_path",
    "is_tool_load_error",
    "load_tool_sources_from_path",
    "yield_tool_sources",
    "yield_tool_sources_on_paths",
]
