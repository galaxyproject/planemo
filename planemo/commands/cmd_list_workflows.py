"""Module describing the planemo ``list_workflows`` command."""

import json

import click

from planemo import options
from planemo.cli import command_function
from planemo.galaxy import profiles
from planemo.galaxy.api import get_workflows
from planemo.io import (
    error,
    info,
)

try:
    from tabulate import tabulate
except ImportError:
    tabulate = None  # type: ignore


def format_url(url):
    if url == "N/A":
        return ""
    return url


@click.command("list_workflows")
@click.option(
    "--raw",
    is_flag=True,
    help="output will be a json structure.",
    default=False,
)
@options.galaxy_url_option()
@options.galaxy_user_key_option()
@options.profile_option()
@command_function
def cli(ctx, raw, **kwds):
    """
    Display available workflows.
    """
    if not raw:
        info("Looking for workflows...")
    profile = kwds.get("profile")
    if profile is not None:
        profile = profiles.ensure_profile(ctx, profile)
        url = profile["galaxy_url"]
        key = profile["galaxy_admin_key"] or profile["galaxy_user_key"]
    else:
        url = kwds.get("galaxy_url")
        key = kwds.get("galaxy_admin_key") or kwds.get("galaxy_user_key")

    workflows = get_workflows(
        url=url,
        key=key,
    )
    if tabulate is not None:
        if raw:
            print(json.dumps(workflows, indent=4, sort_keys=True))
            return
        else:
            print(
                tabulate(
                    {
                        "Workflow ID": workflows.keys(),
                        "Name": [workflow["name"] for _, workflow in workflows.items()],
                        "Url": [format_url(f"{url}/{workflow['url'].strip('/')}") for _, workflow in workflows.items()],
                        "Repo Url": [format_url(workflow["repo_url"]) for _, workflow in workflows.items()],
                    },
                    headers="keys",
                )
            )
    else:
        error("The tabulate package is not installed, invocations could not be listed correctly.")
        print(json.dumps(workflows, indent=4, sort_keys=True))
    if not raw:
        info(f"{len(workflows)} workflows found.")
    return
