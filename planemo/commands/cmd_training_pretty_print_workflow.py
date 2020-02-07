"""Module describing the planemo ``training_pretty_print_workflow`` command."""

import click
import json

from planemo.cli import command_function


@click.command('training_pretty_print_workflow')
@click.argument(
    "path",
    metavar="PATH",
    type=click.Path(
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
    ),
    default="workflow.ga",
)
@command_function
def cli(ctx, path, **kwds):
    """Pretty-print a workflow in-place."""

    with open(path, 'r') as handle:
        data = json.load(handle)

    with open(path, 'w') as handle:
        json.dump(data, handle, sort_keys=True, indent=2)
