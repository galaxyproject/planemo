"""
"""
import sys

import click

from planemo.cli import pass_context
from planemo import options
from planemo import shed
from planemo.io import error
from planemo.io import info
import json


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

    def create(path):
        repo_id = __find_repository(ctx, tsi, path, **kwds)
        if repo_id is None:
            if __create_repository(ctx, tsi, path, **kwds):
                info("Repository created")
                return 0
            else:
                return 2
        else:
            return 1

    if kwds['recursive']:
        exit_code = shed.for_each_repository(create, path)
    else:
        exit_code = create(path)

    sys.exit(exit_code)


def __find_repository(ctx, tsi, path, **kwds):
    """More advanced error handling for finding a repository by ID
    """
    try:
        kwds = kwds.copy()
        kwds["allow_none"] = True
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
        except Exception:
            error(str(e))
        return None
