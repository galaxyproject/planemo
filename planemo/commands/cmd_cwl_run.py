"""Module describing the planemo ``cwl_run`` command."""
import click
from planemo.cli import command_function
from planemo import options
from planemo import cwl


@click.command('cwl_run')
@options.required_tool_arg()
@options.required_job_arg()
@options.galaxy_serve_options()
@options.galaxy_cwl_root_option()
@options.cwl_conformance_test()
@click.option(
    "--cwl_engine",
    type=click.Choice(["galaxy", "cwltool"]),
    default="galaxy",
    help=("Select an engine to run CWL job using, defaults to Galaxy "
          "but the CWL reference implementation cwltool and be selected "
          "also.")
)
@command_function
def cli(ctx, path, job_path, **kwds):
    """Planemo command for running CWL tools and jobs.

    ::

        % planemo cwl_run cat1-tool.cwl cat-job.json
    """
    engine = kwds.get("cwl_engine", "galaxy")
    if engine == "galaxy":
        return cwl.run_galaxy(ctx, path, job_path, **kwds)
    else:
        return cwl.run_cwltool(ctx, path, job_path, **kwds)
