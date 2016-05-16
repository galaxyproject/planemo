"""Module describing the planemo ``profile_list`` command."""
from __future__ import print_function

import click

from planemo.cli import command_function
from planemo.galaxy import profiles


@click.command('profile_list')
@command_function
def cli(ctx, **kwds):
    """List configured profile names."""
    profile_names = profiles.list_profiles(ctx, **kwds)
    print(profile_names)
