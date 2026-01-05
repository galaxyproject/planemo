"""Module describing the planemo ``workflow_job_init`` command."""

import os

import click
from ruamel.yaml import YAML
from ruamel.yaml.comments import CommentedMap

from planemo import options
from planemo.cli import command_function
from planemo.galaxy.workflows import (
    DOCKSTORE_TRS_BASE,
    get_workflow_from_invocation_id,
    is_trs_identifier,
    job_template_with_metadata,
    new_workflow_associated_path,
    TRS_WORKFLOWS_PREFIX,
)
from planemo.io import can_write_to_path


def _convert_to_commented_map(value):
    """Recursively convert dicts to CommentedMap for proper YAML output."""
    if isinstance(value, dict):
        commented = CommentedMap(value)
        for k, v in value.items():
            if isinstance(v, (dict, list)):
                commented[k] = _convert_to_commented_map(v)
        return commented
    elif isinstance(value, list):
        return [_convert_to_commented_map(item) for item in value]
    return value


def _build_commented_yaml(job, metadata):
    """Build a CommentedMap with metadata comments for each input.

    Uses ruamel.yaml to properly add YAML comments with type, description,
    optionality, default value, format, and collection_type information for each input parameter.
    """
    commented = CommentedMap()

    for label, value in job.items():
        # Convert nested dicts/lists to CommentedMap for proper YAML output
        commented[label] = _convert_to_commented_map(value)

        # Add comment with type, description, optionality, default, format, and collection_type
        meta = metadata.get(label, {})
        input_type = meta.get("type", "")
        input_doc = meta.get("doc", "")
        is_optional = meta.get("optional", False)
        default_value = meta.get("default")
        input_format = meta.get("format", "")
        collection_type = meta.get("collection_type", "")

        if input_type or input_doc or is_optional or default_value is not None or input_format or collection_type:
            comment_parts = []
            if input_type:
                comment_parts.append(f"type: {input_type}")
            if collection_type:
                comment_parts.append(f"collection_type: {collection_type}")
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


def _trs_id_to_job_filename(trs_id: str) -> str:
    """Generate a job filename from a TRS ID.

    Args:
        trs_id: TRS ID in various formats:
            - workflow/github.com/org/repo/name[/version]
            - #workflow/github.com/org/repo/name[/version]
            - trs://workflow/github.com/org/repo/name[/version]
            - https://dockstore.org/api/ga4gh/trs/v2/tools/#workflow/...

    Returns:
        A job filename like "workflow-name-job.yml"
    """
    # Strip common prefixes
    identifier = trs_id
    if identifier.startswith(DOCKSTORE_TRS_BASE):
        identifier = identifier[len(DOCKSTORE_TRS_BASE) :]
    if identifier.startswith(TRS_WORKFLOWS_PREFIX):
        identifier = identifier[len(TRS_WORKFLOWS_PREFIX) :]
    if identifier.startswith("#"):
        identifier = identifier[1:]

    # Parse the TRS ID path: workflow/github.com/org/repo/name[/version]
    parts = identifier.split("/")
    if len(parts) >= 5:
        # Get the workflow name (5th element) and optionally version
        workflow_name = parts[4]
        # Sanitize the name for use as filename
        workflow_name = workflow_name.replace(" ", "-").replace("/", "-")
        return f"{workflow_name}-job.yml"
    else:
        # Fallback to a generic name
        return "workflow-job.yml"


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
        if kwds["from_invocation"]:
            output = new_workflow_associated_path(path_basename, suffix="job")
        elif is_trs_identifier(workflow_identifier):
            # Generate output filename from TRS ID
            # Extract workflow name from TRS ID (e.g., workflow/github.com/org/repo/name -> name-job.yml)
            output = _trs_id_to_job_filename(workflow_identifier)
        else:
            output = new_workflow_associated_path(workflow_identifier, suffix="job")
    if not can_write_to_path(output, **kwds):
        ctx.exit(1)

    commented_job = _build_commented_yaml(job, metadata)
    yaml = YAML()
    yaml.default_flow_style = False
    with open(output, "w") as f_job:
        yaml.dump(commented_job, f_job)
