"""Utilities for dealing with continous integration systems."""

import copy
import glob
import math
import os

import yaml

from planemo import (
    git,
    io,
)
from planemo.shed import REPO_METADATA_FILES
from planemo.tools import (
    is_tool_load_error,
    LOAD_ERROR_MESSAGE,
    yield_tool_sources_on_paths,
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
        extended = kwds.get("extended_git_diff", False)
        if path_type == "repo":
            diff_paths = changed_repos_extended(diff_files) if extended else changed_repos(diff_files)
        else:
            diff_paths = changed_tools_extended(ctx, diff_files, cwd) if extended else diff_files

    unique_paths = {os.path.relpath(p, cwd) for p in raw_paths}
    if diff_paths is not None:
        unique_paths = unique_paths.intersection(diff_paths)
    filtered_paths = sorted(io.filter_paths(unique_paths, cwd=cwd, **filter_kwds))
    excluded_paths = sorted(set(unique_paths) - set(filtered_paths))
    if excluded_paths:
        ctx.log("List of excluded paths: %s" % excluded_paths)

    path_count = len(filtered_paths)
    chunk_size = (1.0 * path_count) / kwds["chunk_count"]
    chunk = kwds["chunk"]

    chunked_paths = []
    for i, path in enumerate(filtered_paths):
        if int(math.floor(i / chunk_size)) == chunk:
            chunked_paths.append(path)

    return chunked_paths


def changed_repos(diff_files):
    """Repositories owning ``diff_files``, found by walking up to a metadata file."""
    diff_paths = set()
    for diff_dir in {os.path.dirname(p) for p in diff_files}:
        diff_path = metadata_file_in_path(diff_dir)
        if diff_path:
            diff_paths.add(diff_path)
    return diff_paths


def changed_repos_extended(diff_files):
    """Repositories at, below, or above the directories holding ``diff_files``.

    Descending is what the plain :func:`changed_repos` cannot do: a file shared
    by several repositories -- a tool collection's macros, say -- lives above
    all of them and owns none.

    Assumes the working directory is the repository root. Descending also means
    a changed file sitting directly in a directory of repositories selects all
    of them, which is why this is opt-in behind ``--extended_git_diff``.
    """
    diff_paths = set()
    for diff_dir in {os.path.dirname(p) for p in diff_files}:
        new_diff_paths = set()
        while diff_dir != "" and len(new_diff_paths) == 0:
            for sub_dir in glob.glob(os.path.join(diff_dir, "**", ""), recursive=True):
                # glob yields directories slash-suffixed; filter_paths intersects
                # against os.path.relpath output, which is not.
                diff_path = metadata_file_in_path(os.path.normpath(sub_dir))
                if diff_path:
                    new_diff_paths.add(diff_path)
            diff_dir = os.path.split(diff_dir)[0]
        diff_paths |= new_diff_paths
    return diff_paths


def changed_tools_extended(ctx, diff_files, cwd):
    """Tools owning ``diff_files``, searching each directory and then its parents.

    Resolves a changed test-data file or helper script to the tool that owns it.
    """
    diff_paths = set()
    scanned = {}
    for diff_dir in {os.path.dirname(p) for p in diff_files}:
        tool_paths = set()
        # The directory the change lives in is searched recursively, so that a
        # file shared by several tools -- a collection's macros -- selects each
        # of them. Ancestors are searched shallowly: there we are looking for
        # the tool that owns a changed test-data file or script, not for every
        # tool that happens to sit somewhere beneath a common parent.
        recursive = True
        while diff_dir != "":
            key = (diff_dir, recursive)
            if key not in scanned:
                scanned[key] = _tools_in_path(ctx, diff_dir, recursive)
            tool_paths, contained_tool_files = scanned[key]
            if contained_tool_files:
                break
            diff_dir = os.path.split(diff_dir)[0]
            recursive = False
        diff_paths |= tool_paths
    return {os.path.relpath(p, cwd) for p in diff_paths}


def _tools_in_path(ctx, path, recursive):
    """Loadable tools in ``path``, and whether any tool file was found there.

    The second element separates "no tools here" from "tools here that do not
    parse". Only the former may send the caller further up the tree: a tool
    broken by the very commit under test would otherwise look like an empty
    directory and escalate the selection to every unrelated sibling.
    """
    if not os.path.isdir(path):
        # git diff reports deleted files, whose directory may be gone.
        return set(), False
    tool_paths = set()
    contained_tool_files = False
    for tool_path, tool_source in yield_tool_sources_on_paths(ctx, [path], recursive=recursive):
        contained_tool_files = True
        if is_tool_load_error(tool_source):
            io.error(LOAD_ERROR_MESSAGE % tool_path)
            continue
        tool_paths.add(tool_path)
    return tool_paths, contained_tool_files


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
