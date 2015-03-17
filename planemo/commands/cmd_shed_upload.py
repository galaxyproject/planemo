"""
"""
import click

from planemo.cli import pass_context
from planemo import options
from planemo import shed
from planemo.io import error
from planemo.io import info
from planemo.io import shell
import json

import fnmatch
import os

tar_path = click.Path(
    exists=True,
    file_okay=True,
    dir_okay=False,
    resolve_path=True,
)


# TODO: Implement alternative tool per repo upload strategy.
# TODO: Use git commit hash and origin to generated commit message.
@click.command("shed_upload")
@options.optional_project_arg(exists=True)
@click.option(
    '--message',
    help="Commit message for tool shed upload."
)
@options.shed_owner_option()
@options.shed_name_option()
@options.shed_target_option()
@click.option(
    '--shed_key',
    help="API key for Tool Shed access (required unless e-mail/pass "
         "specified)."
)
@click.option(
    '--shed_email',
    help="E-mail for Tool Shed auth (required unless shed_key is "
         "specified)."
)
@click.option(
    '--shed_password',
    help="Password for Tool Shed auth (required unless shed_key is "
         "specified)."
)
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
    '-r', '--recursive',
    is_flag=True,
    help="Recursively search for repositories to publish to a tool shed",
)
@pass_context
def cli(ctx, path, **kwds):
    """Handle possible recursion through paths for uploading files to a toolshed
    """
    if kwds['recursive']:
        if kwds['name'] is not None:
            error("--name is incompatible with --recursive")
            return -1
        if kwds['tar'] is not None:
            error("--tar is incompatible with --recursive")
            return -1

        ret_codes = []
        for base_path, dirnames, filenames in os.walk(path):
            for filename in fnmatch.filter(filenames, '.shed.yml'):
                ret_codes.append(
                    __handle_upload(ctx, base_path, **kwds)
                )
        # "Good" returns are Nones, everything else is a -1 and should be
        # passed upwards.
        return None if all(x is None for x in ret_codes) else -1
    else:
        return __handle_upload(ctx, path, **kwds)


def __handle_upload(ctx, path, **kwds):
    """Upload a tool directory as a tarball to a tool shed.
    """
    tar_path = kwds.get("tar", None)
    if not tar_path:
        tar_path = shed.build_tarball(path)
    if kwds["tar_only"]:
        shell("cp %s shed_upload.tar.gz" % tar_path)
        return 0
    tsi = shed.tool_shed_client(ctx, **kwds)
    update_kwds = {}
    message = kwds.get("message", None)
    if message:
        update_kwds["commit_message"] = message
    repo_id = __find_repository(ctx, tsi, path, **kwds)
    if repo_id is None and kwds["force_repository_creation"]:
        repo_id = __create_repository(ctx, tsi, path, **kwds)
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
            error("Could not update %s" % path)
            error(exception_content)
            error(e2.read())
        return -1
    info("Repository %s updated successfully." % path)


def __find_repository(ctx, tsi, path, **kwds):
    """More advanced error handling for finding a repository by ID
    """
    try:
        repo_id = shed.find_repository_id(ctx, tsi, path, **kwds)
        return repo_id
    except Exception as e:
        error("Could not update %s" % path)
        try:
            error(e.read())
        except AttributeError:
            # I've seen a case where the error couldn't be read, so now
            # wrapped in try/except
            error("Could not find repository in toolshed")
    return None


def __create_repository(ctx, tsi, path, **kwds):
    """Wrapper for creating the endpoint if it doesn't exist
    """
    try:
        repo = shed.create_repository(ctx, tsi, path, **kwds)
        return repo['id']
    # Have to catch missing snyopsis/bioblend exceptions
    except Exception as e:
        # TODO: galaxyproject/bioblend#126
        try:
            upstream_error = json.loads(e.read())
            error(upstream_error['err_msg'])
        except Exception as e2:
            error(str(e))
        return None
