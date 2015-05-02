""" Utilities for interacting with git using planemo abstractions.
"""

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
