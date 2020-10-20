"""Module describing the planemo ``cwl_run`` command."""
from __future__ import print_function

import json
import os

import click
from galaxy.util import unicodify

from planemo import options
from planemo.cli import command_function
from planemo.engine import engine_context
from planemo.galaxy.profiles import translate_alias
from planemo.io import warn
from planemo.runnable import for_id, for_path
from planemo.tools import uri_to_path


@click.command('run')
@options.required_runnable_arg()
@options.required_job_arg()
@options.galaxy_run_options()
@options.galaxy_config_options()
@options.enable_cwl_option()
@options.galaxy_cwl_root_option()
@options.run_output_directory_option()
@options.run_output_json_option()
@options.engine_options()
@command_function
def cli(ctx, runnable_identifier, job_path, **kwds):
    """Planemo command for running tools and jobs.

    \b
        % planemo run cat1-tool.cwl cat-job.json
    """
    runnable_identifier = translate_alias(ctx, runnable_identifier, kwds.get('profile'))
    path = uri_to_path(ctx, runnable_identifier)
    if os.path.exists(path):
        runnable = for_path(path)
    else:  # assume galaxy workflow id
        runnable = for_id(runnable_identifier)

    # TODO: do a better test of cwl.
    is_cwl = path.endswith(".cwl")
    kwds["cwl"] = is_cwl
    if kwds.get("engine", None) is None:
        if is_cwl:
            kwds["engine"] = "cwltool"
        elif kwds.get('galaxy_url', None):
            kwds["engine"] = "external_galaxy"
        else:
            kwds["engine"] = "galaxy"
    with engine_context(ctx, **kwds) as engine:
        run_result = engine.run(runnable, job_path)
    if not run_result.was_successful:
        warn("Run failed [%s]" % unicodify(run_result))
        ctx.exit(1)
    outputs_dict = run_result.outputs_dict
    output_json = kwds.get("output_json", None)
    if output_json:
        with open(output_json, "w") as f:
            json.dump(outputs_dict, f)

    return 0
