"""
"""
import json
import sys

import click

from planemo.cli import pass_context
from planemo import options
from planemo import shed
from planemo.io import error
from planemo.io import info
from planemo.io import shell


tar_path = click.Path(
    exists=True,
    file_okay=True,
    dir_okay=False,
    resolve_path=True,
)


@click.command("shed_upload")
@options.optional_project_arg(exists=True)
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
@options.recursive_shed_option()
@pass_context
def cli(ctx, path, **kwds):
    """Handle possible recursion through paths for uploading files to a toolshed
    """
    def upload(realized_repository):
        return __handle_upload(ctx, realized_repository, **kwds)

    exit_code = shed.for_each_repository(upload, path, **kwds)
    sys.exit(exit_code)


def __handle_upload(ctx, realized_repository, **kwds):
    """Upload a tool directory as a tarball to a tool shed.
    """
    path = realized_repository.path
    tar_path = kwds.get("tar", None)
    if not tar_path:
        tar_path = shed.build_tarball(path, **kwds)
    if kwds["tar_only"]:
        suffix = ""
        if realized_repository.multiple:
            name = realized_repository.config["name"]
            suffix = "_%s" % name.replace("-", "_")
        shell("cp %s shed_upload%s.tar.gz" % (tar_path, suffix))
        return 0
    tsi = shed.tool_shed_client(ctx, **kwds)
    update_kwds = {}
    message = kwds.get("message", None)
    if message:
        update_kwds["commit_message"] = message

    # TODO: this needs to use realized repository
    repo_id = realized_repository.find_repository_id(ctx, tsi)
    if repo_id is None and kwds["force_repository_creation"]:
        repo_id = realized_repository.create(ctx, tsi)
    # failing to create the repo, give up
    if repo_id is None:
        return -1
    # TODO: support updating repo information if it changes in the config file

    try:
        tsi.repositories.update_repository(repo_id, tar_path, **update_kwds)
    except Exception as e:
        exception_content = e.read()
        try:
            # Galaxy passes nice JSON messages as their errors, which bioblend
            # blindly returns. Attempt to parse those.
            upstream_error = json.loads(exception_content)
            error(upstream_error['err_msg'])
        except Exception as e2:
            error("Could not update %s" % realized_repository.name)
            error(exception_content)
            error(e2.read())
        return -1
    info("Repository %s updated successfully." % realized_repository.name)
    return 0
