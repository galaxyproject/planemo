"""
"""
import sys

import click

from planemo.cli import pass_context
from planemo import options
from planemo import shed


@click.command("shed_diff")
@options.shed_project_arg()
@options.shed_owner_option()
@options.shed_name_option()
@options.shed_target_option()
@options.shed_fail_fast_option()
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
@options.recursive_shed_option()
@pass_context
def cli(ctx, path, **kwds):
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
    """
    def diff(realized_repository):
        return shed.diff_repo(ctx, realized_repository, **kwds)

    exit_code = shed.for_each_repository(diff, path, **kwds)
    sys.exit(exit_code)
