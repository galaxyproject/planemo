"""Module describing the planemo ``lint`` command."""

import click

from planemo import options
from planemo.cli import (
    command_function,
    PlanemoCliContext,
)
from planemo.tool_lint import (
    build_tool_lint_args,
    lint_tools_on_path,
)


@click.command("lint")
@options.optional_tools_arg(multiple=True, allow_uris=True)
@options.report_level_option()
@options.report_xunit()
@options.fail_level_option()
@options.skip_options()
@options.recursive_option()
@options.lint_planemo_defined_tool_linters_options()
# @click.option(
# "--verify",
# is_flag=True,
# help="If an sha256sum is available, download the entire file AND validate it.",
# default=False,
# )
@command_function
def cli(ctx: PlanemoCliContext, uris, **kwds):
    """Check for common errors and best practices."""
    print("LINT")
    lint_args = build_tool_lint_args(ctx, **kwds)
    exit_code = lint_tools_on_path(ctx, uris, lint_args, recursive=kwds["recursive"])

    # TODO: rearchitect XUnit.
    # if kwds['urls']:
    #         collected_data, url_exit_code = check_urls(ctx, paths, **kwds)
    #         if kwds.get('report_xunit', False):
    #             with open(kwds['report_xunit'], 'w') as handle:
    #                 handle.write(build_report.template_data(
    #                     collected_data, template_name='xunit.tpl'))
    ctx.exit(exit_code)
