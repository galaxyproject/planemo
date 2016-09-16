"""Module describing the planemo ``ci_find_repos`` command."""
from __future__ import print_function

import click

from planemo import options
from planemo.ci import filter_paths, print_path_list
from planemo.cli import command_function
from planemo.shed import find_raw_repositories


@click.command('ci_find_repos')
@options.shed_project_arg()
@options.ci_find_options()
@command_function
def cli(ctx, paths, **kwds):
    """Find all shed repositories in one or more directories.

    Currently, a shed repository is considered a directory with a .shed.yml
    file.
    """
    kwds["recursive"] = True
    repos = find_raw_repositories(ctx, paths, **kwds)
    raw_paths = [r.path for r in repos]
    paths = filter_paths(ctx, raw_paths, path_type="dir", **kwds)
    print_path_list(paths, **kwds)
