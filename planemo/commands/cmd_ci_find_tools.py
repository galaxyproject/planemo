"""Module describing the planemo ``ci_find_tools`` command."""
from __future__ import print_function

import click

from planemo import options
from planemo.ci import filter_paths, print_path_list
from planemo.cli import command_function
from planemo.tools import is_tool_load_error, yield_tool_sources_on_paths


@click.command('ci_find_tools')
@options.shed_project_arg()
@options.ci_find_options()
@command_function
def cli(ctx, paths, **kwds):
    """Find all tools in one or more directories.

    Tools can be chunked up, filtered, etc... to build lists of tools to perform
    operations over for continuous integration operations.
    """
    tool_paths = []
    for (tool_path, tool_source) in yield_tool_sources_on_paths(ctx, paths, recursive=True):
        if is_tool_load_error(tool_source):
            continue
        tool_paths.append(tool_path)

    paths = filter_paths(ctx, tool_paths, path_type="file", **kwds)
    print_path_list(paths, **kwds)
