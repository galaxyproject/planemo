import click

from planemo.cli import pass_context
from planemo import virtualenv

VIRTUALENV_PATH_TYPE = click.Path(
    exists=False,
    writable=True,
    resolve_path=True,
)


@click.command("virtualenv")
@click.argument("virtualenv_path",
                metavar="VIRTUALENV_PATH",
                type=VIRTUALENV_PATH_TYPE)
@pass_context
def cli(ctx, virtualenv_path):
    """Create a virtualenv.

    Use virtualenv as library to create a virtualenv for Galaxy if virtualenv
    is not available on the PATH.
    """
    virtualenv.create_and_exit(virtualenv_path)
