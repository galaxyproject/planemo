"""Module describing the planemo ``output_schema`` command."""

import json

import click

from planemo.output_schemas import (
    SCHEMA_MODELS,
    load_output_schemas,
)


@click.command("output_schema")
@click.option("--format", "output_format", default="json", type=click.Choice(["json"]), show_default=True)
@click.option("--schema", "schema_name", type=click.Choice(sorted(SCHEMA_MODELS)), help="Only export one schema.")
def cli(output_format, schema_name):
    """Export JSON Schemas for Planemo machine-readable outputs."""
    click.echo(json.dumps(load_output_schemas(schema_name), indent=2, sort_keys=True))
