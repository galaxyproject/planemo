from __future__ import absolute_import

from galaxy.tools.lint import lint_tool_source

from planemo.exit_codes import (
    EXIT_CODE_GENERIC_FAILURE,
    EXIT_CODE_OK,
)
from planemo.io import (
    coalesce_return_codes,
    error,
    info,
)
from planemo.tools import (
    is_tool_load_error,
    yield_tool_sources_on_paths,
)

LINTING_TOOL_MESSAGE = "Linting tool %s"


def lint_tools_on_path(ctx, paths, lint_args, **kwds):
    assert_tools = kwds.get("assert_tools", True)
    recursive = kwds.get("recursive", False)
    exit_codes = []
    for (tool_path, tool_xml) in yield_tool_sources_on_paths(ctx, paths, recursive):
        if handle_tool_load_error(tool_path, tool_xml):
            exit_codes.append(EXIT_CODE_GENERIC_FAILURE)
            continue
        info("Linting tool %s" % tool_path)
        if not lint_tool_source(tool_xml, **lint_args):
            error("Failed linting")
            exit_codes.append(EXIT_CODE_GENERIC_FAILURE)
        else:
            exit_codes.append(EXIT_CODE_OK)
    return coalesce_return_codes(exit_codes, assert_at_least_one=assert_tools)


def handle_tool_load_error(tool_path, tool_xml):
    """ Return True if tool_xml is tool load error (invalid XML), and
    print a helpful error message.
    """
    is_error = False
    if is_tool_load_error(tool_xml):
        info("Could not lint %s due to malformed xml." % tool_path)
        is_error = True
    return is_error
