"""Module describing the planemo ``tool_factory`` command."""
import os

import click

from planemo import options
from planemo.cli import command_function
from planemo.galaxy import galaxy_serve
from planemo.runnable import for_paths


@click.command('tool_factory')
@options.optional_tools_arg(multiple=True, allow_uris=True)
@options.galaxy_serve_options()
@options.enable_cwl_option()
@options.galaxy_cwl_root_option()
@command_function
def cli(ctx, **kwds):
    """(Experimental) Launch Galaxy with Tool Factory 2.
    For more information about the Galaxy Tool Factory see the publication
    Creating reusable tools from scripts: the Galaxy Tool Factory by Lazarus
    et. al. (10.1093/bioinformatics/bts573). Available at
    http://www.ncbi.nlm.nih.gov/pubmed/23024011.
    workflow isn't working - not sure how to get it loaded - planemo tries to load it as a tool :(
    """
    mod_dir = os.path.dirname(__file__)
    tf_dir = os.path.join(mod_dir, '..', '..', 'planemo_ext', 'tool_factory_2')
    tf_dir = os.path.abspath(tf_dir)
    runnables = for_paths([tf_dir])
    print(f"## tf_dir = {tf_dir}, runnables={runnables}")
    galaxy_serve(ctx, runnables, **kwds)
