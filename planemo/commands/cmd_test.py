import sys

import click

from planemo.cli import pass_context
from planemo import options
from planemo import galaxy_config

from planemo.galaxy_test import (
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
@pass_context
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
    with galaxy_config.galaxy_config(ctx, paths, **kwds) as config:
        return_value = run_in_config(ctx, config, **kwds)
        if return_value:
            sys.exit(return_value)
