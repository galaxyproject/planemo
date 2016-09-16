"""Module describing the planemo ``tool_factory`` command."""
import os

import click

from planemo import options
from planemo.cli import command_function
from planemo.galaxy import serve
from planemo.runnable import for_paths


@click.command('tool_factory')
@options.galaxy_serve_options()
@command_function
def cli(ctx, **kwds):
    """(Experimental) Launch Galaxy with Tool Factory 2.

    For more information about the Galaxy Tool Factory see the publication
    Creating reusable tools from scripts: the Galaxy Tool Factory by Lazarus
    et. al. (10.1093/bioinformatics/bts573). Available at
    http://www.ncbi.nlm.nih.gov/pubmed/23024011.
    """
    mod_dir = os.path.dirname(__file__)
    tf_dir = os.path.join(mod_dir, '..', '..', 'planemo_ext', 'tool_factory_2')
    runnables = for_paths([tf_dir])
    serve(ctx, runnables, **kwds)
