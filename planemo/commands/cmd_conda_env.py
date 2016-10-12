"""Module describing the planemo ``conda_env`` command."""
from __future__ import print_function

import click

from galaxy.tools.deps import conda_util

from planemo import options
from planemo.cli import command_function
from planemo.conda import build_conda_context, collect_conda_targets
from planemo.io import error
from planemo.io import ps1_for_path


SOURCE_COMMAND = """
PRE_CONDA_PS1=$PS1
source %s %s
if [[ -n $BASH_VERSION ]]; then
    hash -r
elif [[ -n $ZSH_VERSION ]]; then
    rehash
else
echo 'Only bash and zsh are supported'
    return 1
fi
PS1="%s"
echo 'Deactivate environment with conda_env_deactivate'
alias conda_env_deactivate="source %s; %s"
"""


@click.command('conda_env')
@options.optional_tools_arg()
@options.conda_target_options()
# @options.skip_install_option()  # TODO
@command_function
def cli(ctx, path, **kwds):
    """Activate a conda environment for tool.

    Source the output of this command to activate a conda environment for this
    tool.

        % . <(planemo conda_env bowtie2.xml)
        % which bowtie2
        TODO_PLACE_PATH_HERE
    """
    conda_context = build_conda_context(ctx, use_planemo_shell_exec=False, **kwds)
    conda_targets = collect_conda_targets(
        ctx, path, conda_context=conda_context
    )
    installed_conda_targets = conda_util.filter_installed_targets(
        conda_targets, conda_context=conda_context
    )
    env_name, exit_code = conda_util.build_isolated_environment(
        installed_conda_targets, conda_context=conda_context
    )
    if exit_code:
        error("Failed to build environmnt for request.")
        return 1

    ps1 = ps1_for_path(path, base="PRE_CONDA_PS1")
    remove_env = "%s env remove -y --name '%s'" % (
        conda_context.conda_exec, env_name
    )
    deactivate = conda_context.deactivate
    activate = conda_context.activate
    command = SOURCE_COMMAND % (
        activate, env_name, ps1,
        deactivate, remove_env
    )
    print(command)
