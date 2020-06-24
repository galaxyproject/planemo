"""Module describing the planemo ``autoupdate`` command."""
import click

from planemo.exit_codes import (
    EXIT_CODE_GENERIC_FAILURE,
    EXIT_CODE_OK
)
from planemo.io import (
    coalesce_return_codes,
    info
)
from planemo.tools import (
    is_tool_load_error,
    yield_tool_sources_on_paths
)

from planemo import options, autoupdate
from planemo.cli import command_function


@click.command('autoupdate')
@options.optional_tools_arg(multiple=True)
@options.report_level_option()
@options.report_xunit()
@options.fail_level_option()
@options.skip_option()
@options.recursive_option()
@command_function
def cli(ctx, paths, **kwds):
    """Auto-update requirements section if necessary"""
    assert_tools = kwds.get("assert_tools", True)
    recursive = kwds.get("recursive", False)
    exit_codes = []
    for (tool_path, tool_xml) in yield_tool_sources_on_paths(ctx, paths, recursive):
        info("Auto-updating tool %s" % tool_path)
        autoupdate.autoupdate(tool_xml)
        if handle_tool_load_error(tool_path, tool_xml):
            exit_codes.append(EXIT_CODE_GENERIC_FAILURE)
            continue
        else:
            exit_codes.append(EXIT_CODE_OK)
    return coalesce_return_codes(exit_codes, assert_at_least_one=assert_tools)
    ctx.exit()


def handle_tool_load_error(tool_path, tool_xml):
    """ Return True if tool_xml is tool load error (invalid XML), and
    print a helpful error message.
    """
    is_error = False
    if is_tool_load_error(tool_xml):
        info("Could not update %s due to malformed xml." % tool_path)
        is_error = True
    return is_error
