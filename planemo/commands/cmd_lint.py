import sys

import click

from planemo.cli import pass_context
from planemo import options

from planemo.tool_lint import build_lint_args
from planemo.tool_lint import lint_tools_on_path


@click.command('lint')
@options.optional_tools_arg(multiple=True)
@options.report_level_option()
@options.fail_level_option()
@options.skip_option()
@options.lint_xsd_option()
@options.recursive_option()
@pass_context
def cli(ctx, paths, **kwds):
    """Check specified tool(s) for common errors and adherence to best
    practices.
    """
    lint_args = build_lint_args(ctx, **kwds)
    exit = lint_tools_on_path(
        ctx,
        paths,
        lint_args,
        recursive=kwds["recursive"]
    )
    sys.exit(exit)
