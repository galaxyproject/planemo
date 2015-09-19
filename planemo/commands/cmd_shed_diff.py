"""
"""
import sys

import click

from planemo.cli import pass_context
from planemo import options
from planemo import shed
from planemo.reports import build_report


@click.command("shed_diff")
@options.shed_read_options()
@click.option(
    "-o", "--output",
    type=click.Path(file_okay=True, resolve_path=True),
    help="Send diff output to specified file.",
    default=None,
)
@click.option(
    '--shed_target_source',
    help="Source Tool Shed to diff against (will ignore local project info"
         " specified). To compare the main Tool Shed against the test, set"
         " this to testtoolshed.",
    default=None,
)
@click.option(
    "--raw",
    is_flag=True,
    help="Do not attempt smart diff of XML to filter out attributes "
         "populated by the Tool Shed.",
)
@click.option(
    "--report_xunit",
    type=click.Path(file_okay=True, resolve_path=True),
    help="Output diff as a faked XUnit report, useful when you want to "
         "automatically diff repositories and be warned when out-of-date.",
    default=None,
)
@pass_context
def cli(ctx, paths, **kwds):
    """Produce diff between local repository and Tool Shed contents.

    By default, this will produce a diff between this repository and what
    would be uploaded to the Tool Shed with the `shed_upload` command - but
    this command can be made to compare other combinations of repositories.
    Here are some examples::

        % # diff for this repository and the main Tool Shed
        % planemo shed_diff
        % # diff for this repository and the test Tool Shed
        % planemo shed_diff --shed_target testtoolshed
        % # diff for the test Tool Shed and main Tool Shed
        % planemo shed_diff --shed_target_source testtoolshed
        % # diff for two an explicitly specified repositories (ignores
        % # current project's shed YAML file.)
        % planemo shed_diff --owner peterjc --name blast_rbh
            --shed_target_source testtoolshed

    This command will return an exit code of:

    - 0 if there are no detected differences.
    - 1 if there are differences.
    - 2 if the target repository doesn't exist.
    - >200 if there are errors attempting to perform a diff.

    **Warning:** ``shed_diff`` doesn't inspect repository metadata, this
    difference applies only to the file contents of files that would actually be
    uploaded to the repository.
    """

    # In a little bit of cheating, we're defining this variable here to collect
    # a "report" on our shed_diff
    collected_data = {
        'results': {
            'total': 0,
            'errors': 0,
            'failures': 0,
            'skips': 0,
        },
        'tests': [],
    }

    def diff(realized_repository):
        result = shed.diff_repo(ctx, realized_repository, **kwds)
        # Collect data about what happened
        collected_data['results']['total'] += 1
        if result >= 200:
            collected_data['results']['errors'] += 1
        elif result > 0:
            collected_data['results']['failures'] += 1
        collected_data['tests'].append({
            'classname': realized_repository.name,
            'result': result,
        })
        return result

    exit_code = shed.for_each_repository(ctx, diff, paths, **kwds)

    if kwds.get('report_xunit', False):
        with open(kwds['report_xunit'], 'w') as handle:
            handle.write(build_report.template_data(
                collected_data, template_name='diff_xunit.tpl'))

    sys.exit(exit_code)
