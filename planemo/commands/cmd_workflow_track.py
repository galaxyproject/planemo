"""Module describing the planemo ``workflow_track`` command."""

import click

from planemo import options
from planemo.cli import command_function
from planemo.engine.factory import engine_context
from planemo.galaxy.workflow_progress import WorkflowProgress


@click.command("workflow_track")
@options.invocation_target_options()
@command_function
def cli(ctx, invocation_id, **kwds):
    """Run defined tests against existing workflow invocation."""
    with WorkflowProgress() as workflow_progress:
        workflow_progress.add_bars()
        import time

        time.sleep(1)
        new_step = {"state": "new"}
        scheduled_step = {"state": "scheduled"}
        new_steps = [new_step, new_step, new_step]
        one_scheduled_steps = [scheduled_step, new_step, new_step]
        two_scheduled_steps = [scheduled_step, scheduled_step, new_step]
        all_scheduled_steps = [scheduled_step, scheduled_step, scheduled_step]
        state_pairs = [
            ({"state": "new"}, {}),
            ({"state": "ready", "steps": new_steps}, {}),
            ({"state": "ready", "steps": one_scheduled_steps}, {"states": {"new": 1}}),
            ({"state": "ready", "steps": two_scheduled_steps}, {"states": {"new": 2}}),
            ({"state": "ready", "steps": two_scheduled_steps}, {"states": {"new": 1, "running": 1}}),
            ({"state": "ready", "steps": two_scheduled_steps}, {"states": {"new": 1, "ok": 1}}),
            ({"state": "ready", "steps": two_scheduled_steps}, {"states": {"ok": 2}}),
            ({"state": "scheduled", "steps": all_scheduled_steps}, {"states": {"ok": 2, "new": 3}}),
            ({"state": "scheduled", "steps": all_scheduled_steps}, {"states": {"ok": 2, "running": 1, "new": 2}}),
            ({"state": "scheduled", "steps": all_scheduled_steps}, {"states": {"ok": 3, "running": 1, "new": 1}}),
            ({"state": "scheduled", "steps": all_scheduled_steps}, {"states": {"ok": 4, "running": 1}}),
            ({"state": "scheduled", "steps": all_scheduled_steps}, {"states": {"ok": 5}}),
        ]
        for invocation, job_states_summary in state_pairs:
            workflow_progress.handle_invocation(invocation, job_states_summary)
            time.sleep(1)

    with engine_context(ctx, engine="external_galaxy", **kwds) as engine, engine.ensure_runnables_served([]) as config:
        user_gi = config.user_gi
        invocation = user_gi.invocations.show_invocation(invocation_id)
        # https://stackoverflow.com/questions/23113494/double-progress-bar-in-python

    ctx.exit(0)
