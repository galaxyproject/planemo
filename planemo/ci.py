"""Utilities for dealing with continous integration systems."""

from __future__ import print_function

import copy
import glob
import math
import os

import yaml

from planemo import git
from planemo import io
from planemo.shed import REPO_METADATA_FILES
from planemo.tools import (
    is_tool_load_error,
    yield_tool_sources_on_paths
)


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
            diff_paths = changed_repos(diff_files)
        else:
            diff_paths = changed_tools(diff_files, ctx, cwd)

    unique_paths = set(os.path.relpath(p, cwd) for p in raw_paths)
    if diff_paths is not None:
        unique_paths = unique_paths.intersection(diff_paths)
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


def changed_repos(diff_files):
    diff_dirs = set(os.path.dirname(p) for p in diff_files)
    diff_paths = set()
    for diff_dir in diff_dirs:
        new_diff_paths = set()
        while diff_dir != "" and len(new_diff_paths) == 0:
            for sub_dir in glob.glob(diff_dir + "/**/", recursive=True):
                diff_path = metadata_file_in_path(sub_dir)
                if diff_path:
                    new_diff_paths.add(diff_path)
            diff_dir = os.path.split(diff_dir)[0]
        diff_paths |= new_diff_paths
    return diff_paths


def changed_tools(diff_files, ctx, cwd):
    diff_paths = set()
    for diff_file in diff_files:
        diff_dir = os.path.dirname(diff_file)
        # search for tool files in each non-root parent*
        new_diff_paths = set()
        while diff_dir != '' and len(new_diff_paths) == 0:
            for (tool_path, tool_source) in yield_tool_sources_on_paths(ctx, [diff_dir], recursive=True):
                if is_tool_load_error(tool_source):
                    continue
                new_diff_paths.add(tool_path)
            diff_dir = os.path.split(diff_dir)[0]
        diff_paths |= new_diff_paths
    diff_paths = set(os.path.relpath(p, cwd) for p in diff_paths)
    return diff_paths


def metadata_file_in_path(diff_dir):
    while diff_dir:
        for metadata_file in REPO_METADATA_FILES:
            if os.path.isfile(os.path.join(diff_dir, metadata_file)):
                return diff_dir
        diff_dir = os.path.dirname(diff_dir)


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
