"""Module describing the planemo ``workflow_lint`` command."""
import json
import sys

import click

from planemo import options
from planemo.cli import command_function


@click.command('workflow_lint')
@options.required_workflow_arg()
@click.option(
    "--topic_name",
    default=None,
    prompt=False,
    help=("Topic of the tutorial where this workflow belongs to.")
)
@command_function
def cli(ctx, workflow_path, **kwds):
    """Lint workflows based on:\n
        - tags attribute\n
        - annotation attribute\n
        - tools in testtoolshed
    """
    kwds["workflows_from_path"] = True
    topic_name = None
    problems = 0
    output = ["---------------------------------------------------------"]

    if workflow_path.endswith(".ga"):
        with open(workflow_path) as json_file:
            data = json.load(json_file)
            # Checking for 'tags' in workflow if topic is known
            if kwds.get("topic_name") and kwds["topic_name"] not in data['tags']:
                problems += 1
                output.append(
                    "{}. The 'tags' attribute is missing. Please add:".format(str(problems), data['name']))
                output.append('"tags": [' + "\n\t" + '"' + kwds["topic_name"] + '"' + "\n]")

            # Checking for 'tags' in workflow
            elif 'tags' not in data or not data['tags']:
                problems += 1
                output.append(
                    "{}. The 'tags' attribute is missing. Please add:".format(str(problems), data['name']))
                output.append('"tags": [' + "\n\t" + '"' + "<topic_name>" + '"' + "\n]")
                
            # Checking for 'annotation' in workflow
            if 'annotation' not in data or not data['annotation']:
                problems += 1
                output.append(
                    "{}. The 'annotation' attribute is missing. Please add:".format(str(problems)))
                output.append('"annotation": "<title of tutorial>"')
                

            # Checking if there are tools used from the testtoolshed
            for stepnr, step in data['steps'].items():
                if step['tool_id'] and step['type'] == 'tool' and 'testtoolshed.g2.bx.psu.edu' in step['tool_id']:
                    problems += 1
                    output.append("{}. Step {} has a tool from the testtoolshed.".format(str(problems), str(stepnr)))     

            if problems:
                output.insert(1, "Workflow '{}' has {} problems because:".format(data['name'], str(problems)))
                output.append("---------------------------------------------------------\n")
                sys.stderr.write("\n".join(output))
                sys.exit(False)

            print('Workflow is ok')
            sys.exit(True)
