import os

import click

from planemo.cli import pass_context
from planemo import galaxy_config
from planemo import galaxy_run
from planemo import options


@click.command('serve')
@options.optional_tools_arg()
@options.galaxy_root_option()
@options.install_galaxy_option()
@options.test_data_option()
@options.dependency_resolvers_option()
@options.job_config_option()
@options.tool_dependency_dir_option()
@options.brew_dependency_resolution()
@pass_context
def cli(ctx, path, **kwds):
    """Launch a Galaxy instance with the specified tool in the tool panel.

    The Galaxy tool panel will include just the referenced tool or tools (by
    default all the tools in the current working directory) and the upload
    tool.

    planemo will search parent directories to see if any is a Galaxy instance
    - but one can pick the Galaxy instance to use with the ``--galaxy_root``
    option or force planemo to download a disposable instance with the
    ``--install_galaxy`` flag.

    ``planemo`` uses temporarily generated config files and environment
    variables to attempt to shield this execution of Galaxy from manually
    launched runs against that same Galaxy root - but this may not be bullet
    proof yet so please careful and do not try this against production Galaxy
    instances.
    """
    # TODO: Preceate a user.
    # TODO: Setup an admin user.
    # TODO: Pass through more parameters.
    # TODO: Populate test-data directory as FTP directory.
    with galaxy_config.galaxy_config(ctx, path, **kwds) as config:
        # TODO: Allow running dockerized Galaxy here instead.
        run_script = os.path.join(config.galaxy_root, "run.sh")
        server_ini = os.path.join(config.config_directory, "galaxy.ini")
        config.env["GALAXY_CONFIG_FILE"] = server_ini
        cmds = [
            run_script,
        ]
        cmd = "; ".join(cmds)
        action = "Starting galaxy"
        galaxy_run.run_galaxy_command(ctx, cmd, config.env, action)
