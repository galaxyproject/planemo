"""
"""
import sys

import click

from planemo.cli import pass_context
from planemo import options
from planemo import shed


tar_path = click.Path(
    exists=True,
    file_okay=True,
    dir_okay=False,
    resolve_path=True,
)


@click.command("shed_upload")
@options.shed_project_arg()
@click.option(
    '--message',
    help="Commit message for tool shed upload."
)
@options.shed_owner_option()
@options.shed_name_option()
@options.shed_target_option()
@options.shed_key_option()
@options.shed_email_option()
@options.shed_password_option()
@click.option(
    '--tar_only',
    is_flag=True,
    help="Produce tar file for upload but do not publish to a tool shed.",
)
@click.option(
    '--tar',
    help="Specify a pre-existing tar file instead of automatically building "
         "one as part of this command.",
    type=tar_path,
    default=None,
)
@click.option(
    '--force_repository_creation',
    help="If a repository cannot be found for the specified user/repo name "
         "pair, then automatically create the repository in the toolshed.",
    is_flag=True,
    default=False
)
@click.option(
    "--check_diff",
    is_flag=True,
    help="Skip uploading if the shed_diff detects there would be no "
         "'difference' (only attributes populated by the shed would would "
         "be updated.)"
)
@options.recursive_shed_option()
@options.shed_fail_fast_option()
@pass_context
def cli(ctx, path, **kwds):
    """Handle possible recursion through paths for uploading files to a toolshed
    """
    def upload(realized_repository):
        return shed.upload_repository(ctx, realized_repository, **kwds)

    exit_code = shed.for_each_repository(upload, path, **kwds)
    sys.exit(exit_code)
