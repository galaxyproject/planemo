"""Utilities for interacting with git using planemo abstractions."""

import os
import subprocess

from galaxy.util import unicodify

from planemo import io


def git_env_for(path):
    """Setup env dictionary to target specified git repo with git commands."""
    env = os.environ.copy()
    env.update({"GIT_WORK_DIR": path, "GIT_DIR": os.path.join(path, ".git")})
    return env


def ls_remote(ctx, remote_repo):
    """Return a dictionary with refs as key and commits as value."""
    commits_and_refs = io.communicate(
        ["git", "ls-remote", remote_repo],
        stdout=subprocess.PIPE,
    )[0]
    return dict(line.split()[::-1] for line in commits_and_refs.decode("utf-8").splitlines())


def init(ctx, repo_path):
    env = git_env_for(repo_path)
    io.communicate(["git", "init"], env=env)


def add(ctx, repo_path, file_path):
    env = git_env_for(repo_path)
    io.communicate(["git", "add", os.path.relpath(file_path, repo_path)], env=env, cwd=repo_path)


def commit(ctx, repo_path, message=""):
    env = git_env_for(repo_path)
    io.communicate(["git", "commit", "-m", message], env=env)


def push(ctx, repo_path, to=None, branch=None, force=False):
    env = git_env_for(repo_path)
    cmd = ["git", "push"]
    if force:
        cmd += ["--force"]
    if to and branch:
        cmd += ["-u", to, branch]
    io.communicate(cmd, env=env, cwd=repo_path)


def branch(ctx, repo_path, branch, from_branch=None):
    env = git_env_for(repo_path)
    cmd = ["git", "checkout", "-b", branch]
    if from_branch is not None:
        cmd.append(from_branch)
    io.communicate(cmd, env=env)


def checkout(ctx, remote_repo, local_path, branch=None, remote="origin", from_branch="master"):
    """Checkout a new branch from a remote repository."""
    env = git_env_for(local_path)
    if not os.path.exists(local_path):
        io.communicate(command_clone(ctx, remote_repo, local_path))
    else:
        io.communicate(["git", "fetch", remote], env=env)

    if branch:
        io.communicate(["git", "checkout", f"{remote}/{from_branch}", "-b", branch], env=env)
    else:
        io.communicate(["git", "merge", "--ff-only", f"{remote}/{from_branch}"], env=env)


def command_clone(ctx, src, dest, mirror=False, branch=None):
    """Produce a command-line string to clone a repository.

    Take in ``ctx`` to allow more configurability down the road.
    """
    cmd = ["git", "clone"]
    if mirror:
        cmd.append("--mirror")
    if branch is not None:
        cmd.extend(["--branch", branch])
    cmd.extend([src, dest])
    return cmd


def diff(ctx, directory, range):
    """Produce a list of diff-ed files for commit range."""
    cmd_template = "cd '%s' && git diff --name-only '%s' --"
    cmd = cmd_template % (directory, range)
    stdout, _ = io.communicate(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    return [line.strip() for line in unicodify(stdout).splitlines() if line]


def clone(*args, **kwds):
    """Clone a git repository.

    See :func:`command_clone` for description of arguments.
    """
    command = command_clone(*args, **kwds)
    return io.communicate(command)


def rev(ctx, directory):
    """Raw revision for git directory specified.

    Throws ``RuntimeError`` if not a git directory.
    """
    cmd_template = "cd '%s' && git rev-parse HEAD"
    cmd = cmd_template % directory
    stdout, _ = io.communicate(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return unicodify(stdout).strip()


def is_rev_dirty(ctx, directory):
    """Check if specified git repository has uncommitted changes."""
    return io.shell(["git", "diff", "--quiet"], cwd=directory) != 0


def rev_if_git(ctx, directory):
    """Determine git revision (or ``None``)."""
    try:
        the_rev = rev(ctx, directory)
        is_dirty = is_rev_dirty(ctx, directory)
        if is_dirty:
            the_rev += "-dirty"
        return the_rev
    except RuntimeError:
        return None
