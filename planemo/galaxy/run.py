"""Utilities for calling Galaxy scripts."""
import os
import string

from galaxy.tools.deps.commands import shell

from planemo.io import info, shell_join
from planemo.virtualenv import create_command


# Activate galaxy's virtualenv if present (needed for tests say but not for
# server because run.sh does this).
ACTIVATE_COMMAND = "[ -e $GALAXY_VIRTUAL_ENV ] && . $GALAXY_VIRTUAL_ENV/bin/activate"
CREATE_COMMAND_TEMPLATE = string.Template(
    'if [ ! -e $GALAXY_VIRTUAL_ENV ]; then $create_virtualenv; fi',
)
PRINT_VENV_COMMAND = shell_join(
    'echo "Set \$GALAXY_VIRTUAL_ENV to $GALAXY_VIRTUAL_ENV"',
    'if [ -e $GALAXY_VIRTUAL_ENV ]',
    'then echo "Virtual environment directory exists."',
    'else echo "Virtual environment directory does not exist."',
    'fi',
)


# TODO: Mac-y curl variant of this.
DOWNLOAD_GALAXY = (
    "wget https://codeload.github.com/galaxyproject/galaxy/tar.gz/"
)

CACHED_VIRTUAL_ENV_COMMAND = ("if [ -d .venv ] || [ -f dist-eggs.ini ];"
                              " then GALAXY_VIRTUAL_ENV=.venv; "
                              " else GALAXY_VIRTUAL_ENV=%s; fi")
UNCACHED_VIRTUAL_ENV_COMMAND = "GALAXY_VIRTUAL_ENV=.venv"


def setup_venv(ctx, kwds):
    if kwds.get("skip_venv", False):
        return ""

    create_template_params = {
        'create_virtualenv': create_command("$GALAXY_VIRTUAL_ENV")
    }
    return shell_join(
        locate_galaxy_virtualenv(ctx, kwds),
        PRINT_VENV_COMMAND if ctx.verbose else None,
        CREATE_COMMAND_TEMPLATE.safe_substitute(create_template_params),
        PRINT_VENV_COMMAND if ctx.verbose else None,
        ACTIVATE_COMMAND,
    )


def locate_galaxy_virtualenv(ctx, kwds):
    if not kwds.get("no_cache_galaxy", False):
        workspace = ctx.workspace
        shared_venv_path = os.path.join(workspace, "gx_venv")
        venv_command = CACHED_VIRTUAL_ENV_COMMAND % shared_venv_path
    else:
        venv_command = UNCACHED_VIRTUAL_ENV_COMMAND
    return shell_join(
        venv_command,
        "export GALAXY_VIRTUAL_ENV",
    )


def shell_if_wheels(command):
    """ Take a shell command and convert it to shell command that runs
    only if Galaxy is new enough to use wheels.
    """
    return "$(grep -q 'skip-venv' run_tests.sh) && %s" % command


def setup_common_startup_args():
    return _set_variable_if_wheels(
        "COMMON_STARTUP_ARGS", "--dev-wheels"
    )


def _set_variable_if_wheels(var, if_wheels_val, else_val=""):
    var_command = '${var}=${else_val}; '
    var_command += shell_if_wheels(
        '${var}="${if_wheels_val}"; '
    )
    var_command += "export ${var}"
    var_command += '; echo "Set ${var} to ${${var}}"'
    return string.Template(var_command).safe_substitute(
        var=var,
        if_wheels_val=if_wheels_val,
        else_val=else_val,
    )


def run_galaxy_command(ctx, command, env, action):
    """Run Galaxy command with informative verbose logging."""
    message = "%s with command [%s]" % (action, command)
    info(message)
    ctx.vlog("With environment variables:")
    ctx.vlog("============================")
    for key, value in env.items():
        ctx.vlog('%s="%s"' % (key, value))
    ctx.vlog("============================")
    exit_code = shell(command, env=env)
    ctx.vlog("run command exited with return code %s" % exit_code)
    return exit_code

__all__ = [
    "setup_venv",
    "run_galaxy_command",
    "DOWNLOAD_GALAXY",
    "setup_common_startup_args",
]
