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
@pass_context
def cli(ctx, path, **kwds):
    """Create a repository in a Galaxy Tool Shed from a ``.shed.yml`` file.
    """
    tsi = shed.tool_shed_client(ctx, **kwds)

    def create(realized_reposiotry):
        repo_id = realized_reposiotry.find_repository_id(ctx, tsi)
        if repo_id is None:
            if realized_reposiotry.create(ctx, tsi):
                info("Repository created")
                return 0
            else:
                return 2
        else:
            return 1

    exit_code = shed.for_each_repository(ctx, create, path, **kwds)
    sys.exit(exit_code)
