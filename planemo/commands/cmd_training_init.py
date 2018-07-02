"""Module describing the planemo ``training_init`` command."""
import os

import click

from planemo import options
from planemo import training
from planemo.config import planemo_option
from planemo.cli import command_function
from planemo.io import write_file
from planemo.runnable import for_path
from planemo.engine import (
    engine_context,
    is_galaxy_engine,
)



@click.command('training_init')
@options.optional_tools_arg(multiple=True, allow_uris=True)
@options.training_init_options()
# @options.force_option()
# @options.galaxy_serve_options()
@command_function
def cli(ctx, uris, **kwds):
    """Build training template from workflow.
    """
    assert is_galaxy_engine(**kwds)
    kwds["no_dependency_resolution"] = True
    training.init(kwds)

#
    #topic_dir = Path("topics") / Path(args.topic_name)
    #if not topic_dir.is_dir():
    #    print("The topic {} does not exist. It will be created".format(args.topic_name))
    #    create_topic(args, topic_dir, template_dir)
#
    #if args.tutorial_name:
    #    tuto_dir = topic_dir / Path("tutorials") / Path(args.tutorial_name)
    #    if not tuto_dir.is_dir():
    #        template_dir = template_dir / Path("tutorials") / Path("tutorial1")
    #        print("The tutorial {} in topic {} does not exist. It will be created.".format(args.tutorial_name, args.topic_name))
    #        create_tutorial(args, tuto_dir, topic_dir, template_dir)
    #    else:
    #        print("The tutorial {} in topic {} already exists. It will be updated with the other arguments".format(args.tutorial_name, args.topic_name))
    #        update_tutorial(args, tuto_dir, topic_dir)
#
#
#
    #if output is None:
    #    output = os.path.splitext(workflow_path)[0] + ".ga"
#
    #runnable = for_path(workflow_path)
    #with engine_context(ctx, **kwds) as galaxy_engine:
    #    with galaxy_engine.ensure_runnables_served([runnable]) as config:
    #        workflow_id = config.workflow_id(workflow_path)
    #        output_dict = config.gi.workflows.export_workflow_dict(workflow_id)
    #        print(output_dict)
    #        import time
    #        time.sleep(10000)
    #        write_file(output, "Test File Here...", force=force)
#