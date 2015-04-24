"""
"""
import sys

import click

from planemo.cli import pass_context
from planemo import options
from planemo import shed
from planemo.io import info


# TODO: Implement alternative tool per repo upload strategy.
# TODO: Use git commit hash and origin to generated commit message.
@click.command("shed_create")
@options.optional_project_arg(exists=True)
@options.shed_owner_option()
@options.shed_name_option()
@options.shed_target_option()
@options.shed_key_option()
@options.shed_email_option()
@options.shed_password_option()
@options.recursive_shed_option()
@pass_context
def cli(ctx, path, **kwds):
    """Create a repository in the toolshed from a .shed.yml file
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

    exit_code = shed.for_each_repository(create, path, **kwds)
    sys.exit(exit_code)
