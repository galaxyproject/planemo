"""
"""
import socket
import sys

import click

from planemo.cli import pass_context
from planemo import options
from planemo import galaxy_serve
from planemo import shed
from planemo import galaxy_test


@click.command("shed_test")
@options.shed_read_options()
@options.galaxy_target_options()
@options.test_options()
@click.option(
    "--skip_dependencies",
    is_flag=True,
    help="Do not install shed dependencies as part of repository installation."
)
@pass_context
def cli(ctx, paths, **kwds):
    """ Run tests of published shed artifacts.

    This command will start a Galaxy instance configured to target the
    specified shed, find published artifacts (tools and dependencies)
    corresponding to command-line arguments and ``.shed.yml`` file(s),
    install these artifacts, and run the tool tests for these commands.

    This command requires the target to be version 15.07 or newer.
    """
    install_args_list = shed.install_arg_lists(ctx, paths, **kwds)
    port = get_free_port()
    kwds["port"] = port
    return_code = 1
    with galaxy_serve.shed_serve(ctx, install_args_list, **kwds) as config:
        config.kill()
        return_code = galaxy_test.run_in_config(
            ctx,
            config,
            installed=True,
            **kwds
        )
    if return_code:
        sys.exit(return_code)


def get_free_port():
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    sock.bind(('localhost', 0))
    port = sock.getsockname()[1]
    sock.close()
    return port
