import os
import sys

import click

from planemo.cli import pass_context
from planemo.io import info
from planemo import options
from planemo import galaxy_config
from planemo import galaxy_run

from galaxy.tools.deps.commands import shell

RUN_TESTS_CMD = (
    "sh run_tests.sh functional.test_toolbox; "
    "return_code=$?; "
    "cp run_functional_tests.html %s; "
    'sh -c "exit $return_code"'
)


@click.command('test')
@options.optional_tools_arg()
@click.option(
    "--test_output",
    type=click.Path(file_okay=True, resolve_path=True),
    help="Output test report.", default="tool_test_output.html"
)
@click.option(
    "--job_output_files",
    type=click.Path(file_okay=False, resolve_path=True),
    help="Write job outputs to directory.", default=None,
)
@options.galaxy_root_option()
@options.install_galaxy_option()
@options.test_data_option()
@options.dependency_resolvers_option()
@options.job_config_option()
@options.tool_dependency_dir_option()
@options.brew_dependency_resolution()
@pass_context
def cli(ctx, path, **kwds):
    """Run the tests in the specified tool tests in a Galaxy instance.

    All referenced tools (by default all the tools in the current working
    directory) will be tested and the resulted disposited in path specified
    with ``--test_output`` (defaults to tool_test_output.html).

    To run these tests planemo needs a Galaxy instance to utilize, planemo
    will search parent directories to see if any is a Galaxy instance
    - but one can pick the Galaxy instance to use with the --galaxy_root
    option or force planemo to download a disposable instance with the
    ``--install_galaxy`` flag.

    planemo uses temporarily generated config files and environment variables
    to attempt to shield this execution of Galaxy from manually launched runs
    against that same Galaxy root - but this may not be bullet proof yet so
    please careful and do not try this against production Galaxy instances.
    """
    with galaxy_config.galaxy_config(path, for_tests=True, **kwds) as config:
        info("Testing using galaxy_root %s", config.galaxy_root)
        # TODO: Allow running dockerized Galaxy here instead.
        server_ini = os.path.join(config.config_directory, "galaxy.ini")
        config.env["GALAXY_CONFIG_FILE"] = server_ini
        config.env["GALAXY_TEST_VERBOSE_ERRORS"] = "true"
        if kwds["job_output_files"]:
            config.env["GALAXY_TEST_SAVE"] = kwds["job_output_files"]
        cd_to_galaxy_command = "cd %s" % config.galaxy_root
        cmd = "; ".join([
            galaxy_run.DEACTIVATE_COMMAND,
            cd_to_galaxy_command,
            galaxy_run.ACTIVATE_COMMAND,
            RUN_TESTS_CMD % kwds["test_output"],
        ])
        info("Running commands [%s]", cmd)
        if shell(cmd, env=config.env):
            sys.exit(1)
