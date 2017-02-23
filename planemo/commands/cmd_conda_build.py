"""Module describing the planemo ``conda_build`` command."""
from __future__ import print_function

import click

from planemo import options
from planemo.cli import command_function
from planemo.conda import build_conda_context


@click.command('conda_build')
@options.conda_target_options(include_local=False)  # No reason to expose local, we have to use it.
@options.recipe_arg(multiple=True)
@command_function
def cli(ctx, paths, **kwds):
    """Perform conda build with Planemo's conda."""
    # Force conda_use_local for building...
    kwds["conda_use_local"] = True
    conda_context = build_conda_context(ctx, handle_auto_init=True, **kwds)
    build_args = list(paths)
    conda_context.exec_command("build", build_args)
