"""Module describing the planemo ``conda_lint`` command."""
import click

from planemo import options
from planemo import conda_lint
from planemo.cli import command_function


@click.command('conda_lint')
@options.report_level_option()
@options.fail_level_option()
@options.recursive_option(
    "Recursively perform command for nested conda directories.",
)
@options.recipe_arg(multiple=True)
@command_function
def cli(ctx, paths, **kwds):
    """Check conda recipe for common issues.

    Built in large part on the work from the BSD licensed anaconda-verify
    project. For more information on anacoda-verify see:
    https://github.com/ContinuumIO/anaconda-verify.
    """
    exit_code = conda_lint.lint_recipes_on_paths(ctx, paths, **kwds)
    ctx.exit(exit_code)
