"""Module describing the planemo ``rerun`` command."""
from __future__ import print_function

import click

from planemo import options
from planemo.cli import command_function
from planemo.engine import engine_context
from planemo.io import error, info
from planemo.runnable import Rerunnable


@click.command('rerun')
@options.profile_option()
@options.galaxy_url_option()
@options.galaxy_user_key_option()
@click.option(
    "--invocation_id",
    multiple=True,
    help=("Rerun failed jobs associated with this invocation ID.")
)
@click.option(
    "--history_id",
    multiple=True,
    help=("Rerun failed jobs associated with this history ID.")
)
@click.option(
    "--job_id",
    multiple=True,
    help=("Rerun failed jobs specified by this job ID.")
)
@command_function
def cli(ctx, **kwds):
    """Planemo command for rerunning and remapping failed jobs on an external Galaxy server.
    One or more history, invocation or job IDs can be specified, and all associated failed
    jobs will be rerun. Note that when specifying a history ID, a maximum of 1000 errored
    jobs will be rerun.

    \b
        % planemo rerun --invocation_id ... --history_id ... --job_id ... --job_id ...
    """
    if not (kwds.get('invocation_id') or kwds.get('history_id') or kwds.get('job_id')):
        error("Please specify at least one history, invocation or job to rerun.")
        ctx.exit(1)
    # Possible TODO: allow collection IDs to be specified as well
    kwds["engine"] = "external_galaxy"
    rerun_successful = True
    with engine_context(ctx, **kwds) as engine:
        for rerunnable_type, rerunnable_ids in {'invocation': 'invocation_id', 'history': 'history_id', 'job': 'job_id'}.items():
            for rerunnable_id in kwds[rerunnable_ids]:
                rerunnable = Rerunnable(rerunnable_id, rerunnable_type, kwds["galaxy_url"])
                rerun_result = engine.rerun(ctx, rerunnable, **kwds)
                if not rerun_result.was_successful:
                    rerun_successful = False

    if rerun_successful:
        info('All requested jobs were rerun successfully.')
        ctx.exit(0)
    else:
        error('Some of the requested jobs could not be rerun.')
        ctx.exit(1)
