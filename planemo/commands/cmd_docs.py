"""Module describing the planemo ``docs`` command."""
import click

from planemo.cli import pass_context

SYNTAX_URL = "http://planemo.readthedocs.org/en/latest/"


@click.command("syntax")
@pass_context
def cli(ctx, **kwds):
    """Open the Planemo documentation in a web browser."""
    click.launch(SYNTAX_URL)
