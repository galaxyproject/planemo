"""Utilities for calling Galaxy scripts."""
import os
import string

from galaxy.tool_util.deps.commands import shell
from six.moves import shlex_quote

from planemo.io import info, shell_join
from planemo.virtualenv import (
    create_command,
    DEFAULT_PYTHON_VERSION,
)


# Activate galaxy's virtualenv if present (needed for tests say but not for
# server because run.sh does this).
ACTIVATE_COMMAND = (
    'if [ -e "$GALAXY_VIRTUAL_ENV" ]; then . "$GALAXY_VIRTUAL_ENV"/bin/activate; '
    'echo "Activated a virtualenv for Galaxy"; echo "$VIRTUAL_ENV"; '
    'else echo "Failed to activate virtualenv."; fi'
)
CREATE_COMMAND_TEMPLATE = string.Template(
    'if [ ! -e "$GALAXY_VIRTUAL_ENV" ]; then $create_virtualenv; echo "Created virtualenv"; fi',
)
PRINT_VENV_COMMAND = shell_join(
    r'echo "Set \$GALAXY_VIRTUAL_ENV to $GALAXY_VIRTUAL_ENV"',
    ('if [ -e "$GALAXY_VIRTUAL_ENV" ]; '
     'then echo "Virtual environment directory exists."; '
     'else echo "Virtual environment directory does not exist."; fi'),
)


CACHED_VIRTUAL_ENV_COMMAND = ("if [ -d .venv ] || [ -f dist-eggs.ini ]; "
                              "then GALAXY_VIRTUAL_ENV=.venv; "
                              "else GALAXY_VIRTUAL_ENV=%s; fi")
UNCACHED_VIRTUAL_ENV_COMMAND = "GALAXY_VIRTUAL_ENV=.venv"


def setup_venv(ctx, kwds):
    if kwds.get("skip_venv", False):
        return ""

    create_template_params = {
        'create_virtualenv': create_command("$GALAXY_VIRTUAL_ENV", kwds.get('galaxy_python_version'))
    }
    return shell_join(
        locate_galaxy_virtualenv(ctx, kwds),
        PRINT_VENV_COMMAND if ctx.verbose else None,
        CREATE_COMMAND_TEMPLATE.safe_substitute(create_template_params),
        PRINT_VENV_COMMAND if ctx.verbose else None,
        ACTIVATE_COMMAND,
        "bash -c 'pwd'" if ctx.verbose else None,
        "bash -c 'which python'" if ctx.verbose else None,
        "bash -c 'which pip'" if ctx.verbose else None,
        "bash -c 'echo $GALAXY_VIRTUAL_ENV'" if ctx.verbose else None,
        "bash -c 'echo $VIRTUAL_ENV'" if ctx.verbose else None,
        "bash -c 'ls -a'" if ctx.verbose else None,
    )


def locate_galaxy_virtualenv(ctx, kwds):
    if os.environ.get("GALAXY_VIRTUAL_ENV"):
        venv_command = ""
    elif not kwds.get("no_cache_galaxy", False):
        workspace = ctx.workspace
        galaxy_branch = kwds.get("galaxy_branch") or "master"
        shared_venv_path = os.path.join(workspace, "gx_venv")
        galaxy_python_version = kwds.get('galaxy_python_version') or DEFAULT_PYTHON_VERSION
        shared_venv_path = "%s_%s" % (shared_venv_path, galaxy_python_version)
        if galaxy_branch != "master":
            shared_venv_path = "%s_%s" % (shared_venv_path, galaxy_branch)
        venv_command = CACHED_VIRTUAL_ENV_COMMAND % shlex_quote(shared_venv_path)
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
    # info not working in pytest+Github actions the way it did in nose?
    info(message)
    ctx.vlog("With environment variables:")
    ctx.vlog("============================")
    for key, value in env.items():
        ctx.vlog('%s="%s"' % (key, value))
    ctx.vlog("============================")
    exit_code = shell(command, env=env)
    ctx.vlog("run command exited with return code %s" % exit_code)
    return exit_code


__all__ = (
    "setup_venv",
    "run_galaxy_command",
    "setup_common_startup_args",
)
