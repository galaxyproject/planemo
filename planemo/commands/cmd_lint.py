import sys
import traceback
import click

from planemo.cli import pass_context
from planemo.io import info
from planemo.io import error
from planemo import options

from galaxy.tools.loader_directory import load_tool_elements_from_path
from galaxy.tools.lint import lint_xml

SKIP_XML_MESSAGE = "Skipping XML file - does not appear to be a tool %s."
LINTING_TOOL_MESSAGE = "Linting tool %s"


@click.command('lint')
@options.optional_tools_arg()
@click.option(
    '--report_level',
    type=click.Choice(['all', 'warn', 'error']),
    default="all"
)
@click.option(
    '--fail_level',
    type=click.Choice(['warn', 'error']),
    default="warn"
)
@pass_context
def cli(ctx, path, report_level="all", fail_level="warn"):
    """Check specified tool(s) for common errors and adherence to best
    practices.
    """
    exit = 0
    lint_args = dict(level=report_level, fail_level=fail_level)
    tools = load_tool_elements_from_path(path, load_exception_handler)
    for (tool_path, tool_xml) in tools:
        if tool_xml.getroot().tag != "tool":
            if ctx.verbose:
                info(SKIP_XML_MESSAGE % tool_path)
            continue
        info("Linting tool %s" % tool_path)
        if not lint_xml(tool_xml, **lint_args):
            exit = 1
    sys.exit(exit)


def load_exception_handler(path, exc_info):
    error("Error loading tool with path %s" % path)
    traceback.print_exception(*exc_info, limit=1, file=sys.stderr)
