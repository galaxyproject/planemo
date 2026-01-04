"""Module describing the planemo ``workflow_job_init`` command."""

import os

import click
from ruamel.yaml import YAML
from ruamel.yaml.comments import CommentedMap

from planemo import options
from planemo.cli import command_function
from planemo.galaxy.workflows import (
    get_workflow_from_invocation_id,
    job_template_with_metadata,
    new_workflow_associated_path,
)
from planemo.io import can_write_to_path


def _build_commented_yaml(job, metadata):
    """Build a CommentedMap with metadata comments for each input.

    Uses ruamel.yaml to properly add YAML comments with type, description,
    optionality, default value, and format information for each input parameter.
    """
    commented = CommentedMap()

    for label, value in job.items():
        # Convert nested dicts to CommentedMap for proper YAML output
        if isinstance(value, dict):
            commented_value = CommentedMap(value)
            # Handle nested elements list for collections
            if "elements" in value and isinstance(value["elements"], list):
                commented_value["elements"] = [
                    CommentedMap(elem) if isinstance(elem, dict) else elem for elem in value["elements"]
                ]
            commented[label] = commented_value
        else:
            commented[label] = value

        # Add comment with type, description, optionality, default, and format if metadata is available
        meta = metadata.get(label, {})
        input_type = meta.get("type", "")
        input_doc = meta.get("doc", "")
        is_optional = meta.get("optional", False)
        default_value = meta.get("default")
        input_format = meta.get("format", "")

        if input_type or input_doc or is_optional or default_value is not None or input_format:
            comment_parts = []
            if input_type:
                comment_parts.append(f"type: {input_type}")
            if input_format:
                comment_parts.append(f"format: {input_format}")
            if input_doc:
                comment_parts.append(f"doc: {input_doc}")
            if is_optional:
                comment_parts.append("optional: true")
            if default_value is not None:
                comment_parts.append(f"default: {default_value}")
            comment_text = ", ".join(comment_parts)
            commented.yaml_set_comment_before_after_key(label, before=comment_text)

    return commented


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

    commented_job = _build_commented_yaml(job, metadata)
    yaml = YAML()
    yaml.default_flow_style = False
    with open(output, "w") as f_job:
        yaml.dump(commented_job, f_job)
