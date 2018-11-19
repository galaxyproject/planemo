"""Module describing the planemo ``workflow_edit`` command."""
import click

from planemo import options
from planemo.cli import command_function
from planemo.engine import (
    engine_context,
    is_galaxy_engine,
)
from planemo.galaxy.serve import sleep_for_serve
from planemo.runnable import (
    for_path,
)


@click.command('workflow_edit')
@options.required_workflow_arg()
@options.galaxy_serve_options()
@command_function
def cli(ctx, workflow_path, output=None, force=False, **kwds):
    """Open a synchronized Galaxy workflow editor.
    """
    assert is_galaxy_engine(**kwds)

    kwds["workflows_from_path"] = True

    runnable = for_path(workflow_path)
    with engine_context(ctx, **kwds) as galaxy_engine:
        with galaxy_engine.ensure_runnables_served([runnable]) as config:
            workflow_id = config.workflow_id(workflow_path)
            url = "%s/workflow/editor?id=%s" % (config.galaxy_url, workflow_id)
            click.launch(url)
            sleep_for_serve()
