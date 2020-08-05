"""Module describing the planemo ``workflow_job_init`` command."""
import click
import yaml

from planemo import options
from planemo.cli import command_function
from planemo.galaxy.workflows import job_template, new_workflow_associated_path
from planemo.io import can_write_to_path


@click.command('workflow_job_init')
@options.required_workflow_arg()
@options.force_option()
@options.workflow_output_artifact()
@command_function
def cli(ctx, workflow_path, output=None, **kwds):
    """Initialize a Galaxy workflow job description for supplied workflow.

    Be sure to your lint your workflow with ``workflow_lint`` before calling this
    to ensure inputs and outputs comply with best practices that make workflow
    testing easier.

    Jobs can be run with the planemo run command (``planemo run workflow.ga job.yml``).
    Planemo run works with Galaxy tools and CWL artifacts (both tools and workflows)
    as well so this command may be renamed to to job_init at something along those
    lines at some point.
    """
    job = job_template(workflow_path)
    if output is None:
        output = new_workflow_associated_path(workflow_path, suffix="job")
    if not can_write_to_path(output, **kwds):
        ctx.exit(1)
    with open(output, "w") as f_job:
        yaml.dump(job, f_job)
