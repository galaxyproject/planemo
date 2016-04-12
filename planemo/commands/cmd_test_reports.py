"""Module describing the planemo ``test_reports`` command."""
import os

import click

from planemo.cli import command_function
from planemo import io
from planemo import options
from planemo.galaxy.test import StructuredData, handle_reports


@click.command('test_reports')
@options.tool_test_json()
@options.test_report_options()
@command_function
def cli(ctx, path, **kwds):
    """Generate various tool test reports (HTML, text, markdown) from
    structure output from tests (tool_test_output.json).
    """
    if not os.path.exists(path):
        io.error("Failed to tool test json file at %s" % path)
        return 1

    test_data = StructuredData(path)
    handle_reports(ctx, test_data.structured_data, kwds)
