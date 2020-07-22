"""Planemo-specific wrappers around galaxy-tool-util tool functionality."""
from __future__ import absolute_import

import os
import sys
import traceback

from galaxy.tool_util import loader_directory
from galaxy.tool_util.fetcher import ToolLocationFetcher

from planemo.io import error, info

is_tool_load_error = loader_directory.is_tool_load_error
SKIP_XML_MESSAGE = "Skipping XML file - does not appear to be a tool %s."
SHED_FILES = ["tool_dependencies.xml", "repository_dependencies.xml"]
LOAD_ERROR_MESSAGE = "Error loading tool with path %s"


def uri_to_path(ctx, uri):
    """Fetch URI to a local path."""
    fetcher = ToolLocationFetcher()
    return fetcher.to_tool_path(uri)


def uris_to_paths(ctx, uris):
    """Fetch multiple URIs to a local path."""
    fetcher = ToolLocationFetcher()
    paths = []
    for uri in uris:
        path = fetcher.to_tool_path(uri)
        paths.append(path)
    return paths


def yield_tool_sources_on_paths(ctx, paths, recursive=False, yield_load_errors=True, exclude_deprecated=False):
    """Walk paths and yield ToolSource objects discovered."""
    for path in paths:
        for (tool_path, tool_source) in yield_tool_sources(ctx, path, recursive, yield_load_errors):
            if exclude_deprecated and 'deprecated' in tool_path:
                continue
            yield (tool_path, tool_source)


def yield_tool_sources(ctx, path, recursive=False, yield_load_errors=True):
    """Walk single path and yield ToolSource objects discovered."""
    tools = load_tool_sources_from_path(
        path,
        recursive,
        register_load_errors=True,
    )
    for (tool_path, tool_source) in tools:
        if is_tool_load_error(tool_source):
            if yield_load_errors:
                yield (tool_path, tool_source)
            else:
                error(LOAD_ERROR_MESSAGE % tool_path)
            continue

        if not _is_tool_source(ctx, tool_path, tool_source):
            continue
        yield (tool_path, tool_source)


def load_tool_sources_from_path(path, recursive, register_load_errors=False):
    """Generate a list for tool sources found down specified path."""
    return loader_directory.load_tool_sources_from_path(
        path,
        _load_exception_handler,
        recursive=recursive,
        register_load_errors=register_load_errors,
    )


def _load_exception_handler(path, exc_info):
    error(LOAD_ERROR_MESSAGE % path)
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


__all__ = (
    "is_tool_load_error",
    "load_tool_sources_from_path",
    "yield_tool_sources",
    "yield_tool_sources_on_paths",
)
