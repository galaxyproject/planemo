"""Module describing the planemo ``list_invocations`` command."""

import json

import click

from planemo import options
from planemo.cli import command_function
from planemo.galaxy import profiles
from planemo.galaxy.api import (
    get_invocations,
    gi,
)
from planemo.galaxy.workflows import remote_runnable_to_workflow_id
from planemo.io import (
    error,
    info,
)
from planemo.runnable_resolve import for_runnable_identifier

try:
    from tabulate import tabulate
except ImportError:
    tabulate = None  # type: ignore


@click.command("list_invocations")
@click.argument(
    "workflow_identifier",
    type=click.STRING,
    required=False,
    default="",
)
@click.option(
    "--raw",
    is_flag=True,
    help="output will be a json structure.",
    default=False,
)
@click.option(
    "--max_items",
    type=click.INT,
    help="max number of items returned.",
    default=100,
)
@click.option(
    "--offset_items",
    type=click.INT,
    help="skip first X items.",
    default=0,
)
@options.profile_option(required=True)
@command_function
def cli(ctx, workflow_identifier, raw, max_items, offset_items, **kwds):
    """
    Get a list of invocations for a particular workflow ID or alias.
    """
    if not raw:
        info(f"Looking for invocations for workflow {workflow_identifier}...")
    profile = profiles.ensure_profile(ctx, kwds.get("profile"))
    if workflow_identifier:
        runnable = for_runnable_identifier(ctx, workflow_identifier, kwds)
        assert runnable.is_remote_workflow_uri
        workflow_id = remote_runnable_to_workflow_id(runnable)
    else:
        workflow_id = ""
    gi_client = gi(None, profile["galaxy_url"], profile["galaxy_admin_key"] or profile["galaxy_user_key"])
    invocations = get_invocations(
        gi=gi_client,
        workflow_id=workflow_id,
        instance=True,
        max_items=max_items,
        offset_items=offset_items,
    )
    if raw:
        print(json.dumps(invocations, indent=4, sort_keys=True))
        return
    if tabulate is not None:
        state_colors = {
            "ok": "\033[92m",  # green
            "running": "\033[93m",  # yellow
            "error": "\033[91m",  # red
            "paused": "\033[96m",  # cyan
            "deleted": "\033[95m",  # magenta
            "deleting": "\033[95m",  # magenta
            "deleted_new": "\033[95m",  # magenta
            "new": "\033[96m",  # cyan
            "queued": "\033[93m",  # yellow
            "skipped": "\033[90m",  # gray
        }

        grouped_invocations = {}
        workflows = {}
        for inv_id, inv in invocations.items():
            wf_id = inv["workflow_id"]
            if wf_id not in grouped_invocations:
                workflow = gi_client.workflows.show_workflow(workflow_id=wf_id, instance=True)
                workflows[wf_id] = (workflow["name"], workflow["id"])
            grouped_invocations.setdefault(wf_id, {})[inv_id] = inv
        for workflow_id, data in grouped_invocations.items():
            header = f"Workflow: {workflows[workflow_id][0]} : {profile['galaxy_url'].strip('/')}/workflows/run?id={workflows[workflow_id][1]}"
            print(f"\n{header}")
            print(len(header) * "=")
            print(
                tabulate(
                    {
                        "Invocation ID": data.keys(),
                        "Invocation report URL": [
                            "{}/workflows/invocations/report?id={}".format(profile["galaxy_url"].strip("/"), inv_id)
                            for inv_id in data
                        ],
                        "History URL": [
                            "{}/histories/view?id={}".format(
                                profile["galaxy_url"].strip("/"), invocations[inv_id]["history_id"]
                            )
                            for inv_id in data
                        ],
                        "Jobs status": [
                            ", ".join([f"{state_colors[k]}{v} jobs {k}\033[0m" for k, v in inv["states"].items()])
                            for inv in data.values()
                        ],
                    },
                    headers="keys",
                )
            )
    else:
        error("The tabulate package is not installed, invocations could not be listed correctly.")
        print(json.dumps(invocations, indent=4, sort_keys=True))
    if not raw:
        info(f"{len(invocations)} invocations found.")
    return
