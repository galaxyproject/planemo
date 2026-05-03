"""Module describing the planemo ``output_schema`` command."""

import json

import click

from planemo.output_schemas import (
    SCHEMA_MODELS,
    load_output_schemas,
)


@click.command("output_schema")
@click.option("--schema", "schema_name", type=click.Choice(sorted(SCHEMA_MODELS)), help="Only export one schema.")
def cli(schema_name):
    """Export JSON Schemas for Planemo machine-readable outputs.

    Available schemas:

    \b
    cli-command-metadata: output from `planemo cli_metadata --command NAME`.
    cli-metadata: output from `planemo cli_metadata`.
    invocation-download-manifest: output from `planemo invocation_download --output_json PATH`.
    run-outputs: output from `planemo run --output_json PATH`.
    test-report: output from test report JSON writers such as `planemo test --test_output_json PATH`.
    """
    click.echo(json.dumps(load_output_schemas(schema_name), indent=2, sort_keys=True))
