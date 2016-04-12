"""Module describing the planemo ``brew`` command."""
import click

from planemo.cli import command_function
from planemo import options

from galaxy.tools.loader_directory import load_tool_elements_from_path

from galaxy.tools.deps.requirements import parse_requirements_from_xml
from galaxy.tools.deps import brew_exts
from galaxy.tools.deps import brew_util
from galaxy.util import bunch


@click.command('brew')
@options.optional_tools_arg()
@options.brew_option()
@command_function
def cli(ctx, path, brew=None):
    """Install tool requirements using brew. (**Experimental**)

    An experimental approach to versioning brew recipes will be used.
    See full discussion on the homebrew-science issues page here -
    https://github.com/Homebrew/homebrew-science/issues/1191. Information
    on the implementation can be found at
    https://github.com/jmchilton/platform-brew
    until a more permanent project home is setup.
    """
    for (tool_path, tool_xml) in load_tool_elements_from_path(path):
        ctx.log('Brewing requirements from tool %s',
                click.format_filename(tool_path))
        mock_args = bunch.Bunch(brew=brew)
        brew_context = brew_exts.BrewContext(mock_args)
        requirements, containers = parse_requirements_from_xml(tool_xml)

        for recipe_context in brew_util.requirements_to_recipe_contexts(
            requirements,
            brew_context
        ):
            brew_exts.versioned_install(recipe_context)
