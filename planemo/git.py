"""Utilities for interacting with git using planemo abstractions."""
import os
import subprocess

from six import text_type

from planemo import io


def git_env_for(path):
    """Setup env dictionary to target specified git repo with git commands."""
    env = {
        "GIT_WORK_DIR": path,
        "GIT_DIR": os.path.join(path, ".git")
    }
    return env


def checkout(ctx, remote_repo, local_path, branch=None, remote="origin", from_branch="master"):
    """Checkout a new branch from a remote repository."""
    env = git_env_for(local_path)
    if not os.path.exists(local_path):
        io.communicate(command_clone(ctx, remote_repo, local_path))
    else:
        io.communicate(["git", "fetch", remote], env=env)

    if branch:
        io.communicate(["git", "checkout", "%s/%s" % (remote, from_branch), "-b", branch], env=env)
    else:
        io.communicate(["git", "merge", "--ff-only", "%s/%s" % (remote, from_branch)], env=env)


def command_clone(ctx, src, dest, bare=False, branch=None):
    """Produce a command-line string to clone a repository.

    Take in ``ctx`` to allow more configurability down the road.
    """
    bare_arg = ""
    if bare:
        bare_arg = "--bare"
    branch_arg = ""
    if branch is not None:
        branch_arg = "--branch '%s'" % branch
    cmd = "git clone %s %s '%s' '%s'" % (bare_arg, branch_arg, src, dest)
    return cmd


def diff(ctx, directory, range):
    """Produce a list of diff-ed files for commit range."""
    cmd_template = "cd '%s' && git diff --name-only '%s' --"
    cmd = cmd_template % (directory, range)
    stdout, _ = io.communicate(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    return [l.strip() for l in text_type(stdout).splitlines() if l]


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
    stdout, _ = io.communicate(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    return text_type(stdout).strip()


def is_rev_dirty(ctx, directory):
    """Check if specified git repository has uncommitted changes."""
    # TODO: Use ENV instead of cd.
    cmd = "cd '%s' && git diff --quiet" % directory
    return io.shell(cmd) != 0


def rev_if_git(ctx, directory):
    """Determine git revision (or ``None``)."""
    try:
        the_rev = rev(ctx, directory)
        is_dirtry = is_rev_dirty(ctx, directory)
        if is_dirtry:
            the_rev += "-dirty"
        return the_rev
    except RuntimeError:
        return None
