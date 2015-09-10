import os

from planemo.io import info
from planemo.io import error

import planemo.linters.xsd

from planemo.tools import load_tool_elements_from_path, is_tool_load_error
from galaxy.tools.lint import lint_xml

SKIP_XML_MESSAGE = "Skipping XML file - does not appear to be a tool %s."
LINTING_TOOL_MESSAGE = "Linting tool %s"
SHED_FILES = ["tool_dependencies.xml", "repository_dependencies.xml"]


def lint_tools_on_path(ctx, paths, lint_args, **kwds):
    assert_tools = kwds.get("assert_tools", True)
    recursive = kwds.get("recursive", False)
    exit = 0
    valid_tools = 0
    for path in paths:
        for (tool_path, tool_xml) in yield_tool_xmls(ctx, path, recursive):
            if handle_tool_load_error(tool_path, tool_xml):
                exit = 1
                continue
            info("Linting tool %s" % tool_path)
            if not lint_xml(tool_xml, **lint_args):
                error("Failed linting")
                exit = 1
            else:
                valid_tools += 1
    if exit == 0 and valid_tools == 0 and assert_tools:
        exit = 2
    return exit


def handle_tool_load_error(tool_path, tool_xml):
    """ Return True if tool_xml is tool load error (invalid XML), and
    print a helpful error message.
    """
    is_error = False
    if is_tool_load_error(tool_xml):
        info("Could not lint %s due to malformed xml." % tool_path)
        is_error = True
    return is_error


def yield_tool_xmls(ctx, path, recursive=False):
    tools = load_tool_elements_from_path(
        path,
        recursive,
        register_load_errors=True,
    )
    for (tool_path, tool_xml) in tools:
        if is_tool_load_error(tool_xml):
            yield (tool_path, tool_xml)
            continue
        if not _is_tool_xml(ctx, tool_path, tool_xml):
            continue
        yield (tool_path, tool_xml)


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
    xsd = kwds.get("xsd", False)
    if xsd:
        return [planemo.linters.xsd]
    else:
        return []


def _is_tool_xml(ctx, tool_path, tool_xml):
    if os.path.basename(tool_path) in SHED_FILES:
        return False
    if tool_xml.getroot().tag != "tool":
        if ctx.verbose:
            info(SKIP_XML_MESSAGE % tool_path)
        return False
    return True
