"""Planemo-specific wrappers around galaxy-lib tool functionality."""
from __future__ import absolute_import

import sys
import traceback
from planemo.io import error

from galaxy.tools import loader_directory

is_tool_load_error = loader_directory.is_tool_load_error


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


__all__ = [
    "load_tool_elements_from_path",
    "is_tool_load_error",
    "load_tool_sources_from_path",
]
