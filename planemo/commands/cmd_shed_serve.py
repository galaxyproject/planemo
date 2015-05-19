"""
"""
import time
import click

from planemo.cli import pass_context
from planemo import options
from planemo import galaxy_serve
from planemo import shed
from planemo import io


@click.command("shed_serve")
@options.shed_read_options()
@options.galaxy_run_options()
@click.option(
    "--skip_dependencies",
    is_flag=True,
    help="Do not install shed dependencies as part of repository installation."
)
@pass_context
def cli(ctx, paths, **kwds):
    """ Serve a transient Galaxy instance after installing repositories
    from a remote Tool Shed.
    """
    install_args_list = shed.install_arg_lists(ctx, paths, **kwds)
    with galaxy_serve.shed_serve(ctx, install_args_list, **kwds) as config:
        gx_url = "http://localhost:%d/" % config.port
        io.info("Galaxy running with tools installed at %s" % gx_url)
        time.sleep(1000000)
