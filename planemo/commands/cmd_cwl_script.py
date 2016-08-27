"""Module describing the planemo ``cwl_script`` command."""
from __future__ import print_function

import click

from planemo import options
from planemo.cli import command_function
from planemo.cwl import to_script


@click.command("cwl_script")
@options.required_tool_arg()
@options.required_job_arg()
@click.option('--no_container', is_flag=True, default=False)
@click.option('outdir', '--output_dir', type=click.Path())
@click.option('basedir', '--base_dir', type=click.Path(), default=".")
@command_function
def cli(ctx, path, job_path, **kwds):
    """Compile simple CWL workflows to shell script.

    Compiles simple common workflow language (CWL) workflows to a shell
    script.
    """
    to_script(ctx, path, job_path, **kwds)
    return 0
