"""Module describing the planemo ``delete_alias`` command."""
import click

from planemo import options
from planemo.cli import command_function
from planemo.galaxy import profiles
from planemo.io import error, info

try:
    from tabulate import tabulate
except ImportError:
    tabulate = None  # type: ignore


@click.command('delete_alias')
@options.alias_option(required=True)
@options.profile_option(required=True)
@command_function
def cli(ctx, alias, profile, **kwds):
    """
    List aliases for a path or a workflow or dataset ID. Aliases are associated with a particular planemo profile.
    """
    info("Looking for profiles...")
    exit_code = profiles.delete_alias(ctx, alias, profile)
    if exit_code == 0:
        info('Alias {} was successfully deleted from profile {}'.format(alias, profile))
    else:
        error('Alias {} does not exist, so was not deleted from profile {}'.format(alias, profile))

    ctx.exit(exit_code)
    return
