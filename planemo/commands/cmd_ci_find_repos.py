"""Module describing the planemo ``ci_find_repos`` command."""
from __future__ import print_function

import copy
import math
import os

import click

from planemo import options
from planemo.cli import command_function
from planemo import git
from planemo.io import filter_paths, open_file_or_standard_output
from planemo.config import planemo_option
from planemo.shed import find_raw_repositories


@click.command('ci_find_repos')
@options.shed_project_arg()
@planemo_option(
    "--exclude",
    type=click.Path(resolve_path=False),
    multiple=True,
    help="Paths to exclude.",
)
@planemo_option(
    "--exclude_from",
    type=click.Path(exists=True, file_okay=True, dir_okay=False, resolve_path=True),
    multiple=True,
    help="File of paths to exclude.",
)
@planemo_option(
    "--changed_in_commit_range",
    help="Exclude paths unchanged in git commit range.",
)
@planemo_option(
    "--chunk_count",
    type=int,
    help="Split output into chunks of this many item and print --chunk such group.",
    default=1,
)
@planemo_option(
    "--chunk",
    type=int,
    help=("When output is split into --chunk_count groups, output the group 0-indexed"
          "by this option."),
    default=0,
)
@planemo_option(
    "--output",
    help="File to output to, or - for standard output.",
    default="-",
)
@command_function
def cli(ctx, paths, **kwds):
    """Find all shed repositories in one or more directories.

    Currently, a shed repository is considered a directory with a .shed.yml
    file.
    """
    kwds["recursive"] = True
    cwd = os.getcwd()
    repos = find_raw_repositories(ctx, paths, **kwds)

    filter_kwds = copy.deepcopy(kwds)
    changed_in_commit_range = kwds.get("changed_in_commit_range", None)
    diff_dirs = None
    if changed_in_commit_range is not None:
        diff_files = git.diff(ctx, cwd, changed_in_commit_range)
        diff_dirs = sorted(set([os.path.dirname(p) for p in diff_files]))

    unique_paths = sorted(set(map(lambda r: os.path.relpath(r.path, cwd), repos)))
    filtered_paths = filter_paths(unique_paths, cwd=cwd, **filter_kwds)
    if diff_dirs is not None:
        new_filtered_paths = []
        for path in filtered_paths:
            if path in diff_dirs:
                new_filtered_paths.append(path)

        filtered_paths = new_filtered_paths

    path_count = len(filtered_paths)
    chunk_size = ((1.0 * path_count) / kwds["chunk_count"])
    chunk = kwds["chunk"]
    with open_file_or_standard_output(kwds["output"], "w") as f:
        for i, path in enumerate(filtered_paths):
            if int(math.floor(i / chunk_size)) == chunk:
                print(path, file=f)
