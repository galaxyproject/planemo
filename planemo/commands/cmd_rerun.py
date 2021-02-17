"""Module describing the planemo ``rerun`` command."""
from __future__ import print_function

import json

import click
from galaxy.util import unicodify

from planemo import options
from planemo.cli import command_function
from planemo.engine import engine_context
from planemo.galaxy.profiles import translate_alias
from planemo.io import warn
from planemo.runnable import Rerunnable


@click.command('rerun')
@options.required_rerunnable_arg()
@options.profile_option()
@options.galaxy_url_option()
@options.galaxy_user_key_option()
@click.option(
    "--failed",
    is_flag=True,
    help=("Rerun only failed jobs of an invocation, remapping over the previously created datasets.")
)
@click.option(
    "--rerunnable-type",
    default="invocation",
    type=click.Choice([
            "invocation",
            "job"
        ]),
    help=("Whether the rerunnable object is a job or an invocation (default invocation).")
)
@command_function
def cli(ctx, rerunnable_identifier, **kwds):
    """Planemo command for rerunning workflow invocations and jobs on an external Galaxy server.

    \b
        % planemo rerun <invocation_id>
    """
    rerunnable_identifier = translate_alias(ctx, rerunnable_identifier, kwds.get('profile'))
    rerunnable = Rerunnable(rerunnable_identifier, kwds["rerunnable_type"], kwds["galaxy_url"])
    kwds["engine"] = "external_galaxy"
    with engine_context(ctx, **kwds) as engine:
        run_result = engine.rerun(ctx, rerunnable, **kwds)
    if not run_result.was_successful:
        warn("Run failed [%s]" % unicodify(run_result))
        ctx.exit(1)
    outputs_dict = run_result.outputs_dict
    output_json = kwds.get("output_json", None)
    if output_json:
        with open(output_json, "w") as f:
            json.dump(outputs_dict, f)
    return 0
