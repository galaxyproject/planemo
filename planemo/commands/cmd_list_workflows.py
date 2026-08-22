"""Module describing the planemo ``list_workflows`` command."""

import json

import click

from planemo import options
from planemo.cli import command_function
from planemo.galaxy import profiles
from planemo.galaxy.api import get_workflows
from planemo.io import (
    info,
    print_table,
)


def workflow_display_url(galaxy_url, workflow_id, published):
    """Build a browsable URL for a workflow.

    The API URL returned by Galaxy renders JSON rather than something a user can
    usefully open, so point at the editor - or the published view when the
    workflow has been shared.
    """
    base = galaxy_url.rstrip("/")
    if published:
        return f"{base}/published/workflow?id={workflow_id}"
    return f"{base}/workflows/edit?id={workflow_id}"


@click.command("list_workflows")
@click.option(
    "--raw",
    is_flag=True,
    help="output will be a json structure.",
    default=False,
)
@options.galaxy_url_option()
@options.galaxy_admin_key_option()
@options.galaxy_user_key_option()
@options.profile_option()
@command_function
def cli(ctx, raw, **kwds):
    """
    Display available workflows.
    """
    profile = kwds.get("profile")
    if profile is not None:
        profile = profiles.ensure_profile(ctx, profile)
        url = profile["galaxy_url"]
        key = profile["galaxy_admin_key"] or profile["galaxy_user_key"]
    else:
        url = kwds.get("galaxy_url")
        key = kwds.get("galaxy_admin_key") or kwds.get("galaxy_user_key")
    if not url:
        raise click.UsageError(
            "Please specify --galaxy_url or --profile to indicate which Galaxy to list workflows from."
        )

    if not raw:
        info("Looking for workflows...")
    workflows = get_workflows(
        url=url,
        key=key,
    )
    if raw:
        print(json.dumps(workflows, indent=4, sort_keys=True))
        return
    print_table(
        {
            "Workflow ID": list(workflows.keys()),
            "Name": [workflow["name"] for workflow in workflows.values()],
            "Published": ["yes" if workflow["published"] else "no" for workflow in workflows.values()],
            "Url": [
                workflow_display_url(url, workflow_id, workflow["published"])
                for workflow_id, workflow in workflows.items()
            ],
            "Repo Url": [workflow["repo_url"] or "" for workflow in workflows.values()],
        }
    )
    info(f"{len(workflows)} workflows found.")
