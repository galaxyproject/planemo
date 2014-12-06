"""
"""
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
@pass_context
def cli(ctx, path, **kwds):
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
    repo_id = shed.find_repository_id(ctx, tsi, path, **kwds)
    try:
        tsi.repositories.update_repository(repo_id, tar_path, **update_kwds)
    except Exception as e:
        error(e.read())
        return -1
    info("Repository updated successfully.")
