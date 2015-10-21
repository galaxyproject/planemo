import click
from planemo.cli import pass_context
from planemo.io import conditionally_captured_io
from planemo import options
from planemo import galaxy_serve
from planemo import cwl


@click.command('cwl_run')
@options.required_tool_arg()
@options.required_job_arg()
@options.galaxy_serve_options()
@options.galaxy_cwl_root_option()
@options.cwl_conformance_test()
@pass_context
def cli(ctx, path, job_path, **kwds):
    """Planemo command for running CWL tools and jobs.

    ::

        % planemo cwl_run cat1-tool.cwl cat-job.json
    """
    # TODO: serve options aren't exactly right - don't care about
    # port for instance.
    kwds["cwl"] = True
    conformance_test = kwds.get("conformance_test", False)
    with conditionally_captured_io(conformance_test):
        with galaxy_serve.serve_daemon(ctx, [path], **kwds) as config:
            cwl_run = cwl.run_cwl_tool(path, job_path, config, **kwds)
    print(cwl_run.cwl_command_state)
    return 0
