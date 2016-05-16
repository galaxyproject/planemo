"""Module describing the planemo ``profile_create`` command."""
from __future__ import print_function

import click

from planemo.cli import command_function
from planemo import options
from planemo.galaxy import profiles


@click.command('profile_create')
@options.profile_name_argument()
@options.profile_database_options()
@command_function
def cli(ctx, profile_name, **kwds):
    """Create a profile."""
    profiles.create_profile(ctx, profile_name, **kwds)
    print("Profile created.")
