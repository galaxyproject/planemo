"""Module describing the planemo ``shed_test`` command."""
import sys

import click
import requests

from planemo import (
    options,
    shed,
)
from planemo.cli import command_function
from planemo.engine import engine_context
from planemo.galaxy.test.actions import handle_reports_and_summary
from planemo.network_util import get_free_port
from planemo.runnable_resolve import for_runnable_identifier


@click.command("shed_test")
@options.shed_read_options()
@options.galaxy_target_options()
@options.test_options()
@click.option(
    "--skip_dependencies", is_flag=True, help="Do not install shed dependencies as part of repository installation."
)
@command_function
def cli(ctx, paths, **kwds):
    """Run tests of published shed artifacts.

    This command will start a Galaxy instance configured to target the
    specified shed, find published artifacts (tools and dependencies)
    corresponding to command-line arguments and ``.shed.yml`` file(s),
    install these artifacts, and run the tool tests for these commands.

    This command requires the target to be version 15.07 or newer.
    """
    install_args_list = kwds["install_args_list"] = shed.install_arg_lists(ctx, paths, **kwds)
    runnables = install_args_list_to_runnables(ctx, install_args_list, kwds)
    kwds["port"] = get_free_port()
    with engine_context(ctx, **kwds) as engine:
        test_data = engine.test(runnables)
        ctx.vlog("engine.test returning [%s]" % test_data)
        return_code = handle_reports_and_summary(ctx, test_data.structured_data, kwds=kwds)
    sys.exit(return_code)


def install_args_list_to_runnables(ctx, install_args_list, kwds):
    runnables = []
    for repo in install_args_list:
        url = f'{repo["tool_shed_url"]}api/repositories/get_repository_revision_install_info'
        response = requests.get(
            url, params={"name": repo["name"], "owner": repo["owner"], "changeset_revision": repo["changeset_revision"]}
        )
        response.raise_for_status()
        install_info = response.json()
        repository_metadata = install_info[1]
        assert repository_metadata["model_class"] == "RepositoryMetadata"
        for tool in repository_metadata.get("valid_tools", []):
            runnable = for_runnable_identifier(ctx, tool["guid"], kwds)
            runnables.append(runnable)
    return runnables
