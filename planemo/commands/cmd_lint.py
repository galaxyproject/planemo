"""Module describing the planemo ``lint`` command."""
import click

from planemo import options
from planemo.cli import command_function
from planemo.lint import build_lint_args
from planemo.tool_lint import lint_tools_on_path


@click.command('lint')
@options.optional_tools_arg(multiple=True, allow_uris=True)
@options.report_level_option()
@options.report_xunit()
@options.fail_level_option()
@options.skip_option()
@options.lint_xsd_option()
@options.recursive_option()
@click.option(
    "--urls",
    is_flag=True,
    default=False,
    help="Check validity of URLs in XML files",
)
@click.option(
    "--doi",
    is_flag=True,
    default=False,
    help="Check validity of DOIs in XML files",
)
@click.option(
    "--conda_requirements",
    is_flag=True,
    default=False,
    help="Check tool requirements for availability in best practice Conda channels.",
)
@click.option(
    "biocontainer",
    "--biocontainer",
    "--biocontainers",
    is_flag=True,
    default=False,
    help="Check best practice BioContainer namespaces for a container definition applicable for this tool.",
)
# @click.option(
# "--verify",
# is_flag=True,
# help="If an sha256sum is available, download the entire file AND validate it.",
# default=False,
# )
@command_function
def cli(ctx, uris, **kwds):
    """Check for common errors and best practices."""
    lint_args = build_lint_args(ctx, **kwds)
    exit_code = lint_tools_on_path(
        ctx,
        uris,
        lint_args,
        recursive=kwds["recursive"]
    )

    # TODO: rearchitect XUnit.
    # if kwds['urls']:
    #         collected_data, url_exit_code = check_urls(ctx, paths, **kwds)
    #         if kwds.get('report_xunit', False):
    #             with open(kwds['report_xunit'], 'w') as handle:
    #                 handle.write(build_report.template_data(
    #                     collected_data, template_name='xunit.tpl'))
    ctx.exit(exit_code)
