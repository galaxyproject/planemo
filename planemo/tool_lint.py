from __future__ import absolute_import

import os

from planemo.io import info
from planemo.io import error
from planemo.io import coalesce_return_codes
from planemo.exit_codes import (
    EXIT_CODE_OK,
    EXIT_CODE_GENERIC_FAILURE,
)

import planemo.linters.xsd
import planemo.linters.urls

from planemo.tools import (
    load_tool_sources_from_path,
    is_tool_load_error,
)
from galaxy.tools.lint import lint_tool_source

SKIP_XML_MESSAGE = "Skipping XML file - does not appear to be a tool %s."
LINTING_TOOL_MESSAGE = "Linting tool %s"
SHED_FILES = ["tool_dependencies.xml", "repository_dependencies.xml"]


def lint_tools_on_path(ctx, paths, lint_args, **kwds):
    assert_tools = kwds.get("assert_tools", True)
    recursive = kwds.get("recursive", False)
    valid_tools = 0
    exit_codes = []
    for path in paths:
        for (tool_path, tool_xml) in yield_tool_sources(ctx, path, recursive):
            if handle_tool_load_error(tool_path, tool_xml):
                exit_codes.append(EXIT_CODE_GENERIC_FAILURE)
                continue
            info("Linting tool %s" % tool_path)
            if not lint_tool_source(tool_xml, **lint_args):
                error("Failed linting")
                exit_codes.append(EXIT_CODE_GENERIC_FAILURE)
            else:
                valid_tools += 1
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


def build_lint_args(ctx, **kwds):
    report_level = kwds.get("report_level", "all")
    fail_level = kwds.get("fail_level", "warn")
    skip = kwds.get("skip", None)
    if skip is None:
        skip = ctx.global_config.get("lint_skip", "")
        if isinstance(skip, list):
            skip = ",".join(skip)

    skip_types = [s.strip() for s in skip.split(",")]
    lint_args = dict(
        level=report_level,
        fail_level=fail_level,
        extra_modules=_lint_extra_modules(**kwds),
        skip_types=skip_types,
    )
    return lint_args


def _lint_extra_modules(**kwds):
    linters = []
    if kwds.get("xsd", False):
        linters.append(planemo.linters.xsd)

    if kwds.get("urls", False):
        linters.append(planemo.linters.urls)

    return linters


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
