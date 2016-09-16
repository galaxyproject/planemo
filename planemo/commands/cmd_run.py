"""Module describing the planemo ``cwl_run`` command."""
from __future__ import print_function

import json

import click

from planemo import options
from planemo.cli import command_function
from planemo.engine import engine_context
from planemo.io import conditionally_captured_io, warn


@click.command('run')
@options.required_tool_arg()
@options.required_job_arg()
@options.galaxy_run_options()
@options.galaxy_config_options()
@options.galaxy_cwl_root_option()
@options.cwl_conformance_test()
@options.run_output_directory_option()
@options.run_output_json_option()
@options.engine_options()
@command_function
def cli(ctx, path, job_path, **kwds):
    """Planemo command for running tools and jobs.

    ::

        % planemo run cat1-tool.cwl cat-job.json
    """
    kwds["cwl"] = path.endswith(".cwl")
    conformance_test = kwds.get("conformance_test", False)

    with conditionally_captured_io(conformance_test):
        with engine_context(ctx, **kwds) as engine:
            run_result = engine.run(path, job_path)

    if not run_result.was_successful:
        warn("Run failed [%s]" % str(run_result))
        ctx.exit(1)

    if conformance_test:
        if hasattr(run_result, "cwl_command_state"):
            command_state = run_result.cwl_command_state
            dumped_json = json.dumps(command_state)
            if hasattr(run_result, "galaxy_paths"):
                for (local_path, galaxy_path) in run_result.galaxy_paths:
                    dumped_json = dumped_json.replace(galaxy_path, local_path)
            print(dumped_json)
    else:
        outputs_dict = run_result.outputs_dict
        print(outputs_dict)
        output_json = kwds.get("output_json", None)
        if output_json:
            with open(output_json, "w") as f:
                json.dump(outputs_dict, f)

    return 0
