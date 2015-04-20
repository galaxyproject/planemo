import sys

import click

from planemo.cli import pass_context
from planemo import options

from planemo.tool_lint import build_lint_args
from planemo.tool_lint import lint_tools_on_path


@click.command('lint')
@options.optional_tools_arg()
@options.report_level_option()
@options.fail_level_option()
@options.lint_xsd_option()
@pass_context
def cli(ctx, path, **kwds):
    """Check specified tool(s) for common errors and adherence to best
    practices.
    """
    lint_args = build_lint_args(**kwds)
    exit = lint_tools_on_path(ctx, path, lint_args)
    sys.exit(exit)
