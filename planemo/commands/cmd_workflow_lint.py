"""Module describing the planemo ``workflow_lint`` command."""
import click

from planemo import options
from planemo.cli import command_function
from planemo.lint import build_lint_args
from planemo.workflow_lint import lint_workflow_artifacts_on_paths


@click.command('workflow_lint')
@options.optional_tools_or_packages_arg(multiple=True)
@options.report_level_option()
@options.report_xunit()
@options.fail_level_option()
@options.skip_option()
@command_function
def cli(ctx, paths, **kwds):
    """Check workflows for syntax errors and best practices."""
    # Unlike tools, lets just make this recursive by default.
    lint_args = build_lint_args(ctx, **kwds)
    exit_code = lint_workflow_artifacts_on_paths(
        ctx,
        paths,
        lint_args,
    )
    ctx.exit(exit_code)
