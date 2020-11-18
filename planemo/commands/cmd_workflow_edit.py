"""Module describing the planemo ``workflow_edit`` command."""
import os

import click

from planemo import options
from planemo.cli import command_function
from planemo.engine import (
    engine_context,
    is_galaxy_engine,
)
from planemo.galaxy.profiles import translate_alias
from planemo.galaxy.serve import sleep_for_serve
from planemo.runnable import (
    for_id,
    for_path,
)


@click.command('workflow_edit')
@options.required_workflow_arg()
@options.galaxy_serve_options()
@command_function
def cli(ctx, workflow_identifier, output=None, force=False, **kwds):
    """Open a synchronized Galaxy workflow editor.
    """
    assert is_galaxy_engine(**kwds)

    workflow_identifier = translate_alias(ctx, workflow_identifier, kwds.get('profile'))
    if os.path.exists(workflow_identifier):
        runnable = for_path(workflow_identifier)
    else:  # assume galaxy workflow id
        runnable = for_id(workflow_identifier)

    kwds["workflows_from_path"] = True

    with engine_context(ctx, **kwds) as galaxy_engine:
        with galaxy_engine.ensure_runnables_served([runnable]) as config:
            workflow_id = config.workflow_id_for_runnable(runnable)
            url = "%s/workflow/editor?id=%s" % (config.galaxy_url, workflow_id)
            click.launch(url)
            if kwds["engine"] != "external_galaxy":
                sleep_for_serve()
