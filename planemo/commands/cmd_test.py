import os
import sys

import click

from planemo.cli import pass_context
from planemo import options
from planemo import galaxy_config

from planemo.galaxy_test import run_in_config

OUTPUT_DFEAULTS = {
    "output": "tool_test_output.html",
    "output_json": "tool_test_output.json",
    "output_xunit": None,
}


@click.command('test')
@options.optional_tools_arg(multiple=True)
@click.option(
    "--test_output",
    type=click.Path(file_okay=True, resolve_path=True),
    help=("Output test report (HTML - for humans) defaults to "
          "tool_test_output.html."),
    default=None,
)
@click.option(
    "--test_output_xunit",
    type=click.Path(file_okay=True, resolve_path=True),
    help="Output test report (xUnit style - for computers).",
    default=None,
)
@click.option(
    "--test_output_json",
    type=click.Path(file_okay=True, resolve_path=True),
    help=("Output test report (planemo json) defaults to "
          "tool_test_output.json."),
    default=None,
)
@click.option(
    "--job_output_files",
    type=click.Path(file_okay=False, resolve_path=True),
    help="Write job outputs to specified directory.",
    default=None,
)
@click.option(
    "--update_test_data",
    is_flag=True,
    help="Update test-data directory with job outputs (normally written to "
         "directory --job_output_files if specified.)"
)
@click.option(
    "--summary",
    type=click.Choice(['none', 'minimal', 'compact']),
    default="minimal",
    help=("Summary style printed to planemo's standard output (see output "
          "reports for more complete summary). Set to 'none' to disable "
          "completely.")
)
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
@options.shed_dependency_resolution()
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
    for name, default in OUTPUT_DFEAULTS.items():
        _populate_default_output(ctx, name, kwds, default)

    kwds["for_tests"] = True
    with galaxy_config.galaxy_config(ctx, paths, **kwds) as config:
        return_value = run_in_config(ctx, config, **kwds)
        if return_value:
            sys.exit(return_value)


def _populate_default_output(ctx, type, kwds, default):
    kwd_key = "test_%s" % type
    kwd_value = kwds.get(kwd_key, None)
    if kwd_value is None:
        global_config = ctx.global_config
        global_config_key = "default_test_%s" % type
        if global_config_key in global_config:
            default_value = global_config[global_config_key]
        else:
            default_value = default

        if default_value:
            default_value = os.path.abspath(default_value)
        kwds[kwd_key] = default_value
