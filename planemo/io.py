import click
from galaxy.tools.deps import commands


def shell(cmds, **kwds):
    info(cmds)
    return commands.shell(cmds, **kwds)


def info(message, *args):
    if args:
        message = message % args
    click.echo(click.style(message, bold=True, fg='green'))


def error(message, *args):
    if args:
        message = message % args
    click.echo(click.style(message, bold=True, fg='red'), err=True)


def warn(message, *args):
    if args:
        message = message % args
    click.echo(click.style(message, fg='red'), err=True)
