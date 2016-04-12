"""Module describing the planemo ``config_init`` command."""
import os
import sys

import click

from planemo.cli import command_function
from planemo import options
from planemo import config
from planemo.io import warn, info

CONFIG_TEMPLATE = """## Planemo Global Configuration File.
## Everything in this file is completely optional - these values can all be
## configured via command line options for the corresponding commands.

## Specify a default galaxy_root for test and server commands here.
#galaxy_root: /path/to/galaxy_root
## Username used with toolshed(s).
#shed_username: "<TODO>"
sheds:
  # For each tool shed you wish to target, uncomment key or both email and
  # password.
  toolshed:
    #key: "<TODO>"
    #email: "<TODO>"
    #password: "<TODO>"
  testtoolshed:
    #key: "<TODO>"
    #email: "<TODO>"
    #password: "<TODO>"
  local:
    #key: "<TODO>"
    #email: "<TODO>"
    #password: "<TODO>"
"""
SUCCESS_MESSAGE = (
    "Wrote configuration template to %s, "
    "please open with editor and fill out."
)


@click.command("config_init")
@options.optional_project_arg(exists=None)
@click.option(
    '--template',
    default=None
)
@command_function
def cli(ctx, path, template=None, **kwds):
    """Help initialize global configuration (in home directory) for Planemo.
    """
    # TODO: prompt for values someday.
    config_path = config.global_config_path()
    if os.path.exists(config_path):
        warn("File %s already exists, exiting." % config_path)
        sys.exit(1)
    with open(config_path, "w") as f:
        f.write(CONFIG_TEMPLATE)
        info(SUCCESS_MESSAGE % config_path)
