import sys

import click

from planemo.cli import pass_context
from planemo import options

from planemo.tool_lint import build_lint_args
from planemo.tool_lint import lint_tools_on_path


@click.command('lint')
@options.optional_tools_arg(multiple=True)
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
# @click.option(
# "--verify",
# is_flag=True,
# help="If an sha256sum is available, download the entire file AND validate it.",
# default=False,
# )
@pass_context
def cli(ctx, paths, **kwds):
    """Check specified tool(s) for common errors and adherence to best
    practices.

    With the --urls flag, this command searches for <package>$URL</package>
    and download actions which specify URLs. Each of those are accessed
    individually. By default, this tool requests the first hundred or so bytes of each listed
    URL and validates that a 200 OK was received. In tool XML files, the --urls
    option checks through the help text for mentioned URLs and checks those.
    """
    lint_args = build_lint_args(ctx, **kwds)
    exit = lint_tools_on_path(
        ctx,
        paths,
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
    sys.exit(exit)
