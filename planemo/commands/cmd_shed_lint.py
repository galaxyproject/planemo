import click
import sys

from planemo.cli import pass_context
from planemo import options
from planemo import shed
from planemo import shed_lint


@click.command('shed_lint')
@options.shed_project_arg()
@options.report_level_option()
@options.fail_level_option()
@options.click.option(
    '--tools',
    is_flag=True,
    default=False,
    help=("Lint tools discovered in the process of linting repositories.")
)
@options.lint_xsd_option()
@options.recursive_shed_option()
@options.shed_fail_fast_option()
@pass_context
def cli(ctx, path, **kwds):
    """Check a Tool Shed repository for common problems.
    """
    def lint(realized_repository):
        return shed_lint.lint_repository(ctx, realized_repository, **kwds)

    kwds["fail_on_missing"] = False
    exit_code = shed.for_each_repository(lint, path, **kwds)
    sys.exit(exit_code)
