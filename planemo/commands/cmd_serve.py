"""Module describing the planemo ``serve`` command."""
import click

from planemo.cli import command_function
from planemo.runnable import for_paths
from planemo.galaxy import galaxy_serve
from planemo import options


@click.command('serve')
@options.optional_tools_arg(multiple=True)
@options.galaxy_serve_options()
@options.enable_cwl_option()
@options.galaxy_cwl_root_option()
@command_function
def cli(ctx, paths, **kwds):
    """Launch a Galaxy instance with the specified tool in the tool panel.

    The Galaxy tool panel will include just the referenced tool or tools (by
    default all the tools in the current working directory) and the upload
    tool.

    planemo will search parent directories to see if any is a Galaxy instance
    - but one can pick the Galaxy instance to use with the ``--galaxy_root``
    option or force planemo to download a disposable instance with the
    ``--install_galaxy`` flag.

    ``planemo`` will run the Galaxy instance in an existing virtualenv if one
    exists in a ``.venv`` directory in the specified ``--galaxy_root``.
    Otherwise, the Galaxy instance will run in a clean virtualenv created in
    ``/tmp``.

    ``planemo`` uses temporarily generated config files and environment
    variables to attempt to shield this execution of Galaxy from manually
    launched runs against that same Galaxy root - but this may not be bullet
    proof yet so please careful and do not try this against production Galaxy
    instances.
    """
    runnables = for_paths(paths)
    galaxy_serve(ctx, runnables, **kwds)
