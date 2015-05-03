""" Utilities for interacting with git using planemo abstractions.
"""
import subprocess
from planemo import io


def command_clone(ctx, src, dest, bare=False):
    """ Take in ctx to allow more configurability down the road.
    """
    bare_arg = ""
    if bare:
        bare_arg = "--bare"
    return "git clone %s '%s' '%s'" % (bare_arg, src, dest)


def clone(*args, **kwds):
    command = command_clone(*args, **kwds)
    return io.shell(command)


def rev(ctx, directory):
    """ Raw revision for git directory specified.

    Throws ``RuntimeError`` if not a git directory.
    """
    cmd_template = "cd '%s' && git rev-parse HEAD"
    cmd = cmd_template % directory
    stdout, _ = io.communicate(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    return stdout.strip()


def is_rev_dirty(ctx, directory):
    cmd = "cd '%s' && git diff --quiet" % directory
    return io.shell(cmd) != 0


def rev_if_git(ctx, directory):
    try:
        the_rev = rev(ctx, directory)
        is_dirtry = is_rev_dirty(ctx, directory)
        if is_dirtry:
            the_rev += "-dirty"
        return the_rev
    except RuntimeError:
        return None
