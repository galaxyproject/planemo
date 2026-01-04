"""Module describing the planemo ``workflow_job_init`` command."""

import os

import click
import yaml

from planemo import options
from planemo.cli import command_function
from planemo.galaxy.workflows import (
    get_workflow_from_invocation_id,
    job_template_with_metadata,
    new_workflow_associated_path,
)
from planemo.io import can_write_to_path


def _dump_with_comments(job, metadata):
    """Dump job template as YAML with metadata comments for each input.

    Since PyYAML doesn't support comments, we manually construct the YAML
    output with comment lines for type and description.
    """
    lines = []
    for label, value in job.items():
        # Add comment with type and description if metadata is available
        meta = metadata.get(label, {})
        input_type = meta.get("type", "")
        input_doc = meta.get("doc", "")

        if input_type or input_doc:
            comment_parts = []
            if input_type:
                comment_parts.append(f"type: {input_type}")
            if input_doc:
                comment_parts.append(f"doc: {input_doc}")
            lines.append(f"# {', '.join(comment_parts)}")

        # Serialize this single key-value pair
        single_item = {label: value}
        yaml_str = yaml.dump(single_item, default_flow_style=False)
        lines.append(yaml_str.rstrip())

    return "\n".join(lines) + "\n"


@click.command("workflow_job_init")
@options.required_workflow_arg()
@options.force_option()
@options.workflow_output_artifact()
@options.galaxy_url_option()
@options.galaxy_user_key_option()
@options.from_invocation()
@options.profile_option()
@command_function
def cli(ctx, workflow_identifier, output=None, **kwds):
    """Initialize a Galaxy workflow job description for supplied workflow.

    Be sure to your lint your workflow with ``workflow_lint`` before calling this
    to ensure inputs and outputs comply with best practices that make workflow
    testing easier.

    Jobs can be run with the planemo run command (``planemo run workflow.ga job.yml``).
    Planemo run works with Galaxy tools and CWL artifacts (both tools and workflows)
    as well so this command may be renamed to to job_init at something along those
    lines at some point.
    """
    if kwds["from_invocation"]:
        if not os.path.isdir("test-data"):
            ctx.log("Creating test-data directory.")
            os.makedirs("test-data")
        path_basename = get_workflow_from_invocation_id(
            workflow_identifier, kwds["galaxy_url"], kwds["galaxy_user_key"]
        )

    job, metadata = job_template_with_metadata(workflow_identifier, **kwds)

    if output is None:
        output = new_workflow_associated_path(
            path_basename if kwds["from_invocation"] else workflow_identifier, suffix="job"
        )
    if not can_write_to_path(output, **kwds):
        ctx.exit(1)
    with open(output, "w") as f_job:
        f_job.write(_dump_with_comments(job, metadata))
