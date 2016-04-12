"""Module describing the planemo ``brew_env`` command."""
from __future__ import print_function
import click
import os

from planemo.cli import command_function
from planemo import options
from planemo.io import ps1_for_path

from galaxy.tools.loader import load_tool
from galaxy.tools.deps.requirements import parse_requirements_from_xml
from galaxy.tools.deps import brew_exts
from galaxy.tools.deps import brew_util
from galaxy.util import bunch


@click.command('brew_env')
@options.optional_tools_arg()
@options.brew_option()
@options.skip_install_option()
@click.option(
    "--shell",
    is_flag=True
)
@command_function
def cli(ctx, path, brew=None, skip_install=False, shell=None):
    """Display commands used to modify environment to inject tool's brew
    dependencies.::

        % . <(planemo brew_env bowtie2.xml)
        % which bowtie2
        /home/john/.linuxbrew/Cellar/bowtie2/2.1.0/bin/bowtie2

    By default this will attempt to attempt to install these recipes as needed.
    This automatic installation can be skipped with the ``--skip_install``
    flag.

    Intead of injecting the enviornment into your current shell using the above
    idiom, the ``--shell`` flag can be sent to launch a new subshell when
    sourced.::

        % . <(planemo brew_env --skip_install --shell bowtie2.xml)
        (bowtie2) % which bowtie2
        /home/john/.linuxbrew/Cellar/bowtie2/2.1.0/bin/bowtie2

    """
    tool_xml = load_tool(path)
    mock_args = bunch.Bunch(brew=brew)
    brew_context = brew_exts.BrewContext(mock_args)
    requirements, containers = parse_requirements_from_xml(tool_xml)

    lines = []
    for recipe_context in brew_util.requirements_to_recipe_contexts(
        requirements,
        brew_context
    ):
        if not skip_install:
            brew_exts.versioned_install(recipe_context)
        lines = brew_exts.build_env_statements_from_recipe_context(
            recipe_context
        )
        split_lines = lines.split("\n")
        lines.extend(split_lines)
    if shell:
        # TODO: Would be cool if this wasn't a bunch of random hackery.
        launch_shell = os.environ.get("SHELL")
        if "bash" in launch_shell:
            ps1 = ps1_for_path(path)
            launch_shell = '(source ~/.bashrc; env PS1="%s" %s --norc)' % (
                ps1,
                launch_shell,
            )
        lines.extend([launch_shell])
        print(";".join(lines))
    else:
        print("\n".join(lines))
