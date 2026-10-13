"""Module describing the planemo ``test_interactive`` command."""

import json

import click

from planemo import options
from planemo.cli import command_function
from planemo.engine.test import test_interactive_runnables
from planemo.galaxy.test import (
    is_interactive_tool,
    test_case_count,
    tool_id_for_path,
)
from planemo.runnable import (
    for_paths,
    RunnableType,
)


@click.command("test_interactive")
@options.optional_tools_arg(multiple=True)
@options.galaxy_target_options()
@options.galaxy_config_options()
@options.run_engine_option()
@options.interactive_test_options()
@command_function
def cli(ctx, paths, **kwds):
    """Test that interactive tools start up and open their entry-point ports.

    An interactive tool, when run, launches a container that opens one or more
    network entry points (ports). This command launches each referenced
    interactive tool on a managed Galaxy instance and verifies that its entry
    points become active and are reachable through the Galaxy Interactive
    Tools (GxIT) proxy.

    By default the test passes if the tool starts and its entry point responds.
    Use ``--expect_content`` to additionally assert that the content served by
    the (first) entry point matches an expected value exactly.

    Interactive tools that require parameters must be given those inputs to
    start at all. Use ``--test_index`` to run a tool with the inputs of one of
    its embedded ``<tests>`` cases (counting starts at 1), or ``--inputs`` to
    supply inputs as a JSON object (which override the test case's inputs).

    As with :func:`planemo.commands.cmd_serve.cli <serve>`, planemo will search
    parent directories to find a Galaxy instance - or force planemo to download
    a disposable one with the ``--install_galaxy`` flag. Interactive tools run
    in containers, so a container runtime (e.g. Docker) must be available.
    """
    runnables = for_paths(paths)
    interactive_runnables = [
        runnable
        for runnable in runnables
        if runnable.type == RunnableType.galaxy_tool and runnable.has_path and is_interactive_tool(runnable.path)
    ]
    skipped = len(runnables) - len(interactive_runnables)
    if skipped:
        ctx.log(f"Skipping {skipped} non-interactive tool(s); only interactive tools are tested.")
    if not interactive_runnables:
        raise click.ClickException("No interactive tools found in the supplied path(s).")

    # Validate --test_index against each tool's embedded test cases up front so
    # we fail fast (with a clear message) before starting a managed Galaxy.
    test_index = kwds.get("test_index")
    if test_index is not None:
        for runnable in interactive_runnables:
            count = test_case_count(runnable.path)
            if not 1 <= test_index <= count:
                raise click.ClickException(
                    f"Tool [{tool_id_for_path(runnable.path)}] has {count} test case(s); "
                    f"--test_index {test_index} is out of range."
                )

    # Interactive tools can only be tested on a managed Galaxy instance.
    engine = kwds.get("engine") or "galaxy"
    if engine not in ("galaxy", "docker_galaxy"):
        raise click.ClickException(f"Engine [{engine}] cannot run interactive tools; use 'galaxy' or 'docker_galaxy'.")
    kwds["engine"] = engine

    # Parse the optional --inputs JSON object.
    inputs = None
    if kwds.get("inputs"):
        try:
            inputs = json.loads(kwds["inputs"])
        except json.JSONDecodeError as e:
            raise click.ClickException(f"Invalid --inputs JSON: {e}")
        if not isinstance(inputs, dict):
            raise click.ClickException('--inputs must be a JSON object (e.g. \'{"param": "value"}\').')
    kwds["inputs"] = inputs

    return_value = test_interactive_runnables(ctx, interactive_runnables, **kwds)
    ctx.exit(return_value)
