"""
"""
import click

from planemo.cli import pass_context

SYNTAX_URL = "https://wiki.galaxyproject.org/Admin/Tools/ToolConfigSyntax"


@click.command("syntax")
@pass_context
def cli(ctx, **kwds):
    """Open tool config syntax wiki page in a web browser.
    """
    click.launch(SYNTAX_URL)
