"""Module describing the planemo ``docs`` command."""
import click

from planemo.cli import command_function

SYNTAX_URL = "http://planemo.readthedocs.org/en/latest/"


@click.command("syntax")
@command_function
def cli(ctx, **kwds):
    """Open the Planemo documentation in a web browser."""
    click.launch(SYNTAX_URL)
