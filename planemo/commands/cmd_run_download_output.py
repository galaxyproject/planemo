"""Module describing the planemo ``run`` command."""

import click

from bioblend.galaxy import GalaxyInstance
from bioblend.galaxy.invocations import InvocationClient

from planemo import options
from planemo.galaxy.activity import invocation_to_run_response
from planemo.cli import command_function
from planemo.galaxy import profiles
from planemo.io import info
from planemo.runnable_resolve import for_runnable_identifier
import os


@click.command("run_download_output")
@options.run_output_directory_option()
@options.run_output_json_option()
@options.profile_option(required=True)
@click.option(
    "-i",
    "--invocation",
    help="invocation identifier",
    type=click.STRING,
)
@click.option(
    "-w",
    "--workflow",
    help="workflow identifier",
    type=click.STRING,
)
@click.option(
    "-h",
    "--history",
    help="history identifier",
    type=click.STRING,
)
@click.option(
    "--ignore_missing_output",
    is_flag=True,
    default=True,
    help="Ignore missing output files"
)
@command_function
def cli(ctx, invocation, workflow, history, output_directory, ignore_missing_output, **kwds):
    """Download output files from a completed Galaxy workflow invocation.

    This command allows downloading output files after a workflow has been executed
    through Galaxy's web interface or through planemo.

    \b
        % planemo run_download_output --profile dev --invocation <invocation_id> --workflow <workflow_id> --history <history_id>
    """

    profile = profiles.ensure_profile(ctx, kwds.get("profile"))
    gi = GalaxyInstance(url=profile["galaxy_url"], key=profile["galaxy_admin_key"] or profile["galaxy_user_key"])
    invocation_client = InvocationClient(gi)

    invocation_data = invocation_client.show_invocation(invocation)
    # NEED I way of decoding workflow id and history_id from the invocation data,
    # we have the information but need to decode it before requesting the files.
    invocation_data['workflow_id'] = workflow
    invocation_data['history_id'] = history
    runnable = for_runnable_identifier(ctx, invocation_data['workflow_id'], kwds)

    run_response = invocation_to_run_response(ctx, gi, runnable, invocation_data)

    if output_directory is None:
        output_directory = f"output_{invocation}"
        os.makedirs(output_directory, exist_ok=True)
        info(f"No output_directory provided, saving to {output_directory}")
    run_response.collect_outputs(output_directory, ignore_missing_output=ignore_missing_output)
