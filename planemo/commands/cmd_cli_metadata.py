"""Module describing the planemo ``cli_metadata`` command."""

import json

import click

from planemo.cli_metadata import (
    load_command_metadata,
    load_planemo_metadata,
)


@click.command("cli_metadata")
@click.option("--format", "output_format", default="json", type=click.Choice(["json"]), show_default=True)
@click.option("--command", "command_name", help="Only export metadata for the selected command.")
@click.option("--include-internal", is_flag=True, help="Include internal commands not documented as public API.")
def cli(output_format, command_name, include_internal):
    """Export structured metadata for Planemo CLI commands."""
    if command_name:
        metadata = load_command_metadata(command_name)
    else:
        metadata = load_planemo_metadata(include_internal=include_internal)
    click.echo(json.dumps(metadata, indent=2, sort_keys=True))
