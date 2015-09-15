"""
"""
import sys

import click

from planemo.cli import pass_context
from planemo import options
from planemo import shed
from planemo.io import info, error


@click.command("shed_update")
@options.shed_publish_options()
@options.shed_upload_options()
@options.shed_skip_upload()
@options.shed_skip_metadata()
@pass_context
def cli(ctx, paths, **kwds):
    """Update repository in shed from a ``.shed.yml`` file.

    By default this command will update both repository metadata
    from ``.shed.yml`` and upload new contents from the repository
    directory.

    ::

        % planemo shed_update

    This will update the main tool shed with the repository defined
    by a ``.shed.yml`` file in the current working directory. Both
    the location of the ``.shed.yml`` and the tool shed to upload to
    can be easily configured. For instance, the following command can
    be used if ``.shed.yml`` if contained in ``path/to/repo`` and the
    desire is to update the test tool shed.

    ::

        % planemo shed_update --shed_target testtoolshed path/to/repo

    Another important option is ``--check_diff`` - this doesn't affect the
    updating of shed metadata but it will check for differences before
    uploading new contents to the tool shed. This may important because the
    tool shed will automatically populate certain attributes in tool shed
    artifact files (such as ``tool_dependencies.xml``) and this may
    cause unwanted installable revisions to be created when there are no
    important changes.

    The lower-level ``shed_upload`` command should be used instead if
    the repository doesn't define complete metadata in a ``.shed.yml``.
    """
    tsi = shed.tool_shed_client(ctx, **kwds)

    def update(realized_repository):
        upload_ret_code = 0
        upload_ok = True
        if not kwds["skip_upload"]:
            upload_ret_code = shed.upload_repository(
                ctx, realized_repository, **kwds
            )
            upload_ok = not upload_ret_code
        if upload_ret_code == 2:
            error("Failed to update repository it does not exist "
                  "in target ToolShed.")
            return upload_ret_code
        repo_id = realized_repository.find_repository_id(ctx, tsi)
        metadata_ok = True
        if not kwds["skip_metadata"]:
            metadata_ok = realized_repository.update(ctx, tsi, repo_id)
        if metadata_ok:
            info("Repository metadata updated.")
        else:
            error("Failed to update repository metadata.")
        if metadata_ok and upload_ok:
            return 0
        else:
            error("Failed to update a repository.")
            return 1

    exit_code = shed.for_each_repository(ctx, update, paths, **kwds)
    sys.exit(exit_code)
