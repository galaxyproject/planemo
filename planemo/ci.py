"""Utilities for dealing with continous integration systems."""

from __future__ import print_function

import copy
import math
import os

import yaml

from planemo import git
from planemo import io
from planemo.shed import SHED_CONFIG_NAME


def filter_paths(ctx, raw_paths, path_type="repo", **kwds):
    """Filter ``paths``.

    ``path_type`` is ``repo`` or ``file``.
    """
    cwd = os.getcwd()

    filter_kwds = copy.deepcopy(kwds)
    changed_in_commit_range = kwds.get("changed_in_commit_range", None)
    diff_paths = None
    if changed_in_commit_range is not None:
        diff_files = git.diff(ctx, cwd, changed_in_commit_range)
        if path_type == "repo":
            diff_dirs = set(os.path.dirname(p) for p in diff_files)
            diff_paths = set()
            for diff_dir in diff_dirs:
                while diff_dir:
                    if os.path.isfile(os.path.join(diff_dir, SHED_CONFIG_NAME)):
                        diff_paths.add(diff_dir)
                        break
                    diff_dir = os.path.dirname(diff_dir)
        else:
            diff_paths = diff_files

    unique_paths = set(os.path.relpath(p, cwd) for p in raw_paths)
    if diff_paths is not None:
        new_unique_paths = []
        for path in unique_paths:
            if path in diff_paths:
                new_unique_paths.append(path)
        unique_paths = new_unique_paths
    filtered_paths = sorted(io.filter_paths(unique_paths, cwd=cwd, **filter_kwds))
    excluded_paths = sorted(set(unique_paths) - set(filtered_paths))
    if excluded_paths:
        ctx.log("List of excluded paths: %s" % excluded_paths)

    path_count = len(filtered_paths)
    chunk_size = ((1.0 * path_count) / kwds["chunk_count"])
    chunk = kwds["chunk"]

    chunked_paths = []
    for i, path in enumerate(filtered_paths):
        if int(math.floor(i / chunk_size)) == chunk:
            chunked_paths.append(path)

    return chunked_paths


def group_paths(paths):
    repos = {}
    for path in paths:
        repo = os.path.split(path)[0]
        if repo not in repos:
            repos[repo] = []
        repos[repo].append(path)
    return [" ".join(repos[_]) for _ in repos]


def print_path_list(paths, **kwds):
    with io.open_file_or_standard_output(kwds["output"], "w") as f:
        for path in paths:
            print(path, file=f)


def print_as_yaml(item, **kwds):
    with io.open_file_or_standard_output(kwds["output"], "w") as f:
        f.write(yaml.safe_dump(item))
