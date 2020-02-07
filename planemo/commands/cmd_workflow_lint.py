"""Module describing the planemo ``workflow_lint`` command."""
import json
import os
import sys
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


@click.command('workflow_lint')
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
    """Lint workflows.
    """
    assert is_galaxy_engine(**kwds)

    kwds["no_dependency_resolution"] = True

    if workflow_path.endswith(".ga"):
        import sys 
        import json 

        with open(workflow_path) as json_file:
            data = json.load(json_file)
            # Checking for 'tags' in workflow
            if 'tags' not in data:
                sys.stderr.write(
                    "Workflow {} has no corresponding 'tags' attribute.")
                sys.exit(False)

            # Checking for 'annotation' in workflow
            elif 'annotation' not in data or not data['annotation']:
                sys.stderr.write(
                    "Workflow {} has no corresponding 'annotation' attribute. Please add: \n".format(data['name']))
                sys.stderr.write('"annotation": "<title of tutorial>"' + "\n")
                sys.exit(False)

            # Checking if there are tools used from the testtoolshed
            else:
                for stepnr, step in data['steps'].items():
                    if step['tool_id'] and step['type'] == 'tool' and 'testtoolshed.g2.bx.psu.edu' in step['tool_id']:
                        sys.stderr.write("Workflow {} has a tool from the testtoolshed in step {}.\n".format(
                            data['name'], str(stepnr)))
                        sys.exit(False)
                print('Workflow is ok')

        