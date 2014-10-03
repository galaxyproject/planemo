import click
from planemo.cli import pass_context


@click.command('tool_init')
@pass_context
def cli(ctx):
    """Not yet implemented, but someone should write a wizard for initializing
    Galaxy tools and open a pull request :).
    """
    raise NotImplementedError()
