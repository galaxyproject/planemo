"""Module describing the planemo ``test`` command."""
import click

from planemo.cli import command_function
from planemo.engine import (
    engine_context,
    is_galaxy_engine,
)
from planemo import options
from planemo.runnable import (
    for_paths,
    RunnableType,
)
from planemo.galaxy import galaxy_config
from planemo.io import info

from planemo.galaxy.test import (
    handle_reports_and_summary,
    run_in_config,
)


@click.command('test')
@options.optional_tools_arg(multiple=True)
@click.option(
    "--failed",
    is_flag=True,
    help="Re-run only failed tests. This command will read "
         "tool_test_output.json to determine which tests failed so this "
         "file must have been produced with the same set of tool ids "
         "previously.",
    default=False,
)
@options.galaxy_target_options()
@options.galaxy_config_options()
@options.test_options()
@options.engine_options()
@command_function
def cli(ctx, paths, **kwds):
    """Run the tests in the specified tool tests in a Galaxy instance.

    All referenced tools (by default all the tools in the current working
    directory) will be tested and the results quickly summarized.

    To run these tests planemo needs a Galaxy instance to utilize, planemo
    will search parent directories to see if any is a Galaxy instance
    - but one can pick the Galaxy instance to use with the --galaxy_root
    option or force planemo to download a disposable instance with the
    ``--install_galaxy`` flag.

    In additon to to quick summary printed to the console - various detailed
    output summaries can be configured. ``tool_test_output.html`` (settable
    via ``--test_output``) will contain a human consumable HTML report
    describing the test run. A JSON file (settable via ``--test_output_json``
    and defaulting to ``tool_test_output.json``) will also be created. These
    files can can be disabled by passing in empty arguments or globally by
    setting the values ``default_test_output`` and/or
    ``default_test_output_json`` in ``~/.planemo.yml`` to ``null``. For
    continuous integration testing a xUnit-style report can be confiured using
    the ``--test_output_xunit``.

    planemo uses temporarily generated config files and environment variables
    to attempt to shield this execution of Galaxy from manually launched runs
    against that same Galaxy root - but this may not be bullet proof yet so
    please careful and do not try this against production Galaxy instances.
    """
    kwds["for_tests"] = True
    runnables = for_paths(paths)
    enable_beta_test = any([r.type not in [RunnableType.galaxy_tool, RunnableType.directory] for r in runnables])
    enable_beta_test = enable_beta_test or not is_galaxy_engine(**kwds)
    if enable_beta_test:
        info("Enable beta testing mode to test artifact that isn't a Galaxy tool.")
        with engine_context(ctx, **kwds) as engine:
            test_data = engine.test(runnables)
            return_value = handle_reports_and_summary(ctx, test_data.structured_data, kwds=kwds)
    else:
        with galaxy_config(ctx, runnables, **kwds) as config:
            return_value = run_in_config(ctx, config, **kwds)

    ctx.exit(return_value)
