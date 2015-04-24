import click
import sys

from planemo.cli import pass_context
from planemo import options
from planemo import shed
from planemo import shed_lint


@click.command('shed_lint')
@options.optional_project_arg(exists=True)
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
@pass_context
def cli(ctx, path, recursive=False, **kwds):
    """Check a Tool Shed repository for common problems.
    """
    def lint(realized_repository):
        path = realized_repository.real_path
        return shed_lint.lint_repository(ctx, path, **kwds)

    exit_code = shed.for_each_repository(lint, path)
    sys.exit(exit_code)
