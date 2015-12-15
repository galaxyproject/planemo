import click

from planemo.cli import pass_context
from planemo.io import shell
from planmo import options
from galaxy.tools.deps import conda_util


@click.command('conda_init')
@options.conda_prefix_option()
@pass_context
def cli(ctx, conda_prefix=None):
    """Download and install conda.

    This will download conda for managing dependencies for your platform
    using the appropriate Miniconda installer.

    By running this command, you are agreeing to the terms of the conda
    license a 3-clause BSD 3 license. Please review full license at
    http://docs.continuum.io/anaconda/eula.
    """
    conda_util(conda_prefix, shell_exec=shell)
