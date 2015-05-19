import click
import sys

from planemo.cli import pass_context
from planemo import options
from planemo import shed
from planemo import shed_lint


@click.command('shed_lint')
@options.shed_realization_options()
@options.report_level_option()
@options.fail_level_option()
@options.click.option(
    '--tools',
    is_flag=True,
    default=False,
    help=("Lint tools discovered in the process of linting repositories.")
)
@options.lint_xsd_option()
@options.click.option(
    '--ensure_metadata',
    is_flag=True,
    default=False,
    help=("Ensure .shed.yml files contain enough metadata for each repository "
          "to allow automated creation and/or updates.")
)
@pass_context
def cli(ctx, paths, **kwds):
    """Check a Tool Shed repository for common problems.
    """
    def lint(realized_repository):
        return shed_lint.lint_repository(ctx, realized_repository, **kwds)

    kwds["fail_on_missing"] = False
    exit_code = shed.for_each_repository(ctx, lint, paths, **kwds)
    sys.exit(exit_code)
