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
    info,
    print_table,
)
from planemo.runnable_resolve import for_runnable_identifier

STATE_COLORS = {
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


def _format_job_state(state, count):
    color = STATE_COLORS.get(state, "")
    reset = "\033[0m" if color else ""
    return f"{color}{count} jobs {state}{reset}"


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
    help="Output invocations as JSON.",
    default=False,
)
@click.option(
    "--max-items",
    "--max_items",
    type=click.IntRange(min=0),
    help="Maximum number of invocations to return.",
    default=100,
    show_default=True,
)
@click.option(
    "--offset-items",
    "--offset_items",
    type=click.IntRange(min=0),
    help="Skip this many invocations before returning results.",
    default=0,
    show_default=True,
)
@options.profile_option(required=True)
@command_function
def cli(ctx, workflow_identifier, raw, max_items, offset_items, **kwds):
    """
    Get invocations, optionally filtering by workflow ID or alias.
    """
    if not raw:
        scope = f"workflow {workflow_identifier}" if workflow_identifier else "all workflows"
        info(f"Looking for invocations for {scope}...")
    profile = profiles.ensure_profile(ctx, kwds.get("profile"))
    if workflow_identifier:
        runnable = for_runnable_identifier(ctx, workflow_identifier, kwds)
        assert runnable.is_remote_workflow_uri
        workflow_id = remote_runnable_to_workflow_id(runnable)
    else:
        workflow_id = None
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

    galaxy_url = profile["galaxy_url"].rstrip("/")
    grouped_invocations = {}
    workflows = {}
    for invocation_id, invocation in invocations.items():
        workflow_id = invocation["workflow_id"]
        if workflow_id not in workflows:
            workflow = gi_client.workflows.show_workflow(workflow_id=workflow_id, instance=True)
            workflows[workflow_id] = workflow
        grouped_invocations.setdefault(workflow_id, {})[invocation_id] = invocation

    for workflow_id, workflow_invocations in grouped_invocations.items():
        workflow = workflows[workflow_id]
        header = f"Workflow: {workflow['name']} : {galaxy_url}/workflows/run?id={workflow['id']}"
        print(f"\n{header}")
        print("=" * len(header))
        print_table(
            {
                "Invocation ID": list(workflow_invocations),
                "Invocation report URL": [
                    f"{galaxy_url}/workflows/invocations/report?id={invocation_id}"
                    for invocation_id in workflow_invocations
                ],
                "History URL": [
                    f"{galaxy_url}/histories/view?id={invocation['history_id']}"
                    for invocation in workflow_invocations.values()
                ],
                "Jobs status": [
                    ", ".join(_format_job_state(state, count) for state, count in invocation["states"].items())
                    for invocation in workflow_invocations.values()
                ],
            }
        )
    info(f"{len(invocations)} invocations found.")
    return
