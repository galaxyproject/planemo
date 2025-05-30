"""Module describing the planemo ``download_run_output`` command."""

import os

import click
from bioblend.galaxy import GalaxyInstance

from planemo import options
from planemo.cli import command_function
from planemo.galaxy import profiles
from planemo.galaxy.activity import invocation_to_run_response
from planemo.io import info
from planemo.runnable_resolve import for_runnable_identifier


@click.command("download_run_output")
@options.run_output_directory_option()
@options.run_output_json_option()
@options.profile_option(required=True)
@options.galaxy_url_option()
@options.galaxy_user_key_option()
@click.argument("invocation_id", type=click.STRING)
@click.option("--ignore_missing_output", is_flag=True, default=True, help="Ignore missing output files")
@command_function
def cli(ctx, invocation, output_directory, ignore_missing_output, **kwds):
    """Download output files from a completed Galaxy workflow invocation.

    This command allows downloading output files after a workflow has been executed
    through Galaxy's web interface or through planemo.

    \b
        % planemo download_run_output INVOCATION_ID
    """

    profile = profiles.ensure_profile(ctx, kwds.get("profile"))
    gi = GalaxyInstance(url=profile["galaxy_url"], key=profile["galaxy_admin_key"] or profile["galaxy_user_key"])

    invocation_data = gi.invocations.show_invocation(invocation)
    workflow_id = gi.workflows.show_workflow(invocation_data["workflow_id"], instance=True)["id"]
    runnable = for_runnable_identifier(ctx, workflow_id, kwds)
    run_response = invocation_to_run_response(ctx, gi, runnable, invocation_data)
    if output_directory is None:
        output_directory = f"output_{invocation}"
        os.makedirs(output_directory, exist_ok=True)
        info(f"No output_directory provided, saving to {output_directory}")
    run_response.collect_outputs(output_directory, ignore_missing_output=ignore_missing_output)
