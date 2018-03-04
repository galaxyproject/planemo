"""Module describing the planemo ``workflow_convert`` command."""
import json
import os

import click


from planemo import options
from planemo.cli import command_function
from planemo.config import planemo_option
from planemo.engine import (
    engine_context,
    is_galaxy_engine,
)
from planemo.io import write_file
from planemo.runnable import (
    for_path,
)


@click.command('workflow_convert')
@options.required_workflow_arg()
@options.force_option()
@planemo_option(
    "-o", "--output",
    default=None,
    type=click.Path(
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
    )
)
@options.galaxy_serve_options()
@command_function
def cli(ctx, workflow_path, output=None, force=False, **kwds):
    """Convert Format 2 workflow to a native Galaxy workflow.
    """
    assert is_galaxy_engine(**kwds)

    kwds["no_dependency_resolution"] = True

    if output is None:
        output = os.path.splitext(workflow_path)[0] + ".ga"

    runnable = for_path(workflow_path)
    with engine_context(ctx, **kwds) as galaxy_engine:
        with galaxy_engine.serve_runnables([runnable]) as config:
            workflow_id = config.workflow_id(workflow_path)
            output_dict = config.gi.workflows.export_workflow_dict(workflow_id)

            output_contents = json.dumps(output_dict)
            write_file(output, output_contents, force=force)
