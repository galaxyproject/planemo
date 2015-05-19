"""
"""
import sys

import click

from planemo.cli import pass_context
from planemo import options
from planemo import shed
from planemo.io import info


@click.command("shed_create")
@options.shed_publish_options()
@options.shed_message_option()
@options.shed_skip_upload()
@pass_context
def cli(ctx, paths, **kwds):
    """Create a repository in a Galaxy Tool Shed from a ``.shed.yml`` file.
    """
    tsi = shed.tool_shed_client(ctx, **kwds)

    def create(realized_repository):
        repo_id = realized_repository.find_repository_id(ctx, tsi)
        if repo_id is None:
            if realized_repository.create(ctx, tsi):
                info("Repository created")
                if not kwds["skip_upload"]:
                    return shed.upload_repository(
                        ctx, realized_repository, **kwds
                    )
                else:
                    return 0
            else:
                return 2
        else:
            return 1

    exit_code = shed.for_each_repository(ctx, create, paths, **kwds)
    sys.exit(exit_code)
