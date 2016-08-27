"""Module describing the planemo ``brew_init`` command."""
import urllib
from tempfile import mkstemp

import click

from planemo.cli import command_function
from planemo.io import shell


INSTALL_SCRIPT = "https://raw.github.com/Homebrew/linuxbrew/go/install"


@click.command('brew_init')
@command_function
def cli(ctx):
    """Download linuxbrew install & run it in ruby.

    Linuxbrew is a fork of Homebrew (http://brew.sh/linuxbrew/).

    For more information on installing linuxbrew and pre-requisites see
    https://github.com/Homebrew/linuxbrew#installation.

    Homebrew or linuxbrew are required in order to use the other commands
    ``brew`` and ``brew_shell``.
    """
    fname = mkstemp('install_brew')
    urllib.urlretrieve(INSTALL_SCRIPT, fname)
    shell(["ruby", fname])
