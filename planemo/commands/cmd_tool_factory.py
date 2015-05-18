import os

import click
from planemo.cli import pass_context
from planemo import options
from planemo import galaxy_serve


@click.command('tool_factory')
@options.galaxy_serve_options()
@pass_context
def cli(ctx, **kwds):
    """(Experimental) Launch Galaxy with the Tool Factory 2 available.

    For more information about the Galaxy Tool Factory see the publication
    Creating reusable tools from scripts: the Galaxy Tool Factory by Lazarus
    et. al. (10.1093/bioinformatics/bts573). Available at
    http://www.ncbi.nlm.nih.gov/pubmed/23024011.
    """
    mod_dir = os.path.dirname(__file__)
    tf_dir = os.path.join(mod_dir, '..', '..', 'planemo_ext', 'tool_factory_2')
    galaxy_serve.serve(ctx, [os.path.abspath(tf_dir)], **kwds)
