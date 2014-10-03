import click

import urllib
from tempfile import mktemp

from planemo.cli import pass_context
from galaxy.tools.deps import commands


INSTALL_SCRIPT = "https://raw.github.com/Homebrew/linuxbrew/go/install"


@click.command('brew_init')
@pass_context
def cli(ctx):
    """Download linuxbrew install and run it with ruby. Linuxbrew is a fork
    of Homebrew (http://brew.sh/linuxbrew/).

    For more information on installing linuxbrew and pre-requisites see
    https://github.com/Homebrew/linuxbrew#installation.

    Homebrew or linuxbrew are required in order to use the other commands
    ``brew`` and ``brew_shell``.
    """
    fname = mktemp('install_brew')
    urllib.urlretrieve(INSTALL_SCRIPT, fname)
    commands.execute(["ruby", fname])
