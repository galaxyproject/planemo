"""Module describing the planemo ``shed_serve`` command."""
import time

import click

from planemo import io
from planemo import options
from planemo import shed
from planemo.cli import command_function
from planemo.galaxy import shed_serve


@click.command("shed_serve")
@options.shed_read_options()
@options.galaxy_run_options()
@options.galaxy_config_options()
@options.pid_file_option()
@click.option(
    "--skip_dependencies",
    is_flag=True,
    help="Do not install shed dependencies as part of repository installation."
)
@command_function
def cli(ctx, paths, **kwds):
    """Launch Galaxy with Tool Shed dependencies.

    This command will start a Galaxy instance configured to target the
    specified shed, find published artifacts (tools and dependencies)
    corresponding to command-line arguments and ``.shed.yml`` file(s),
    install these artifacts, and serve a Galaxy instances that can be
    logged into and explored interactively.
    """
    install_args_list = shed.install_arg_lists(ctx, paths, **kwds)
    with shed_serve(ctx, install_args_list, **kwds) as config:
        gx_url = "http://localhost:%d/" % config.port
        io.info("Galaxy running with tools installed at %s" % gx_url)
        time.sleep(1000000)
