"""Module describing the planemo ``cwl_run`` command."""
from __future__ import print_function

import json

import click
from galaxy.util import unicodify

from planemo import options
from planemo.cli import command_function
from planemo.engine import engine_context
from planemo.io import warn
from planemo.runnable_resolve import for_runnable_identifier


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
    runnable = for_runnable_identifier(ctx, runnable_identifier, kwds)

    is_cwl = runnable.type.is_cwl_artifact
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
