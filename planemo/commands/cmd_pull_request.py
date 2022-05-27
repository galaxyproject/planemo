"""Module describing the planemo ``pull_request`` command."""
import click

from planemo import (
    github_util,
    options,
)
from planemo.cli import command_function
from planemo.config import planemo_option


@click.command("pull_request")
@planemo_option(
    "-m", "--message", type=click.STRING, default=None, help="Message describing the pull request to create."
)
@options.optional_project_arg(exists=None)
@command_function
def cli(ctx, path, message=None, **kwds):
    """Short-cut to quickly create a pull request for a relevant Github repo.

    For instance, the following will clone, fork, and branch the tools-iuc
    repository to allow a subsequent pull request to fix a problem with bwa.

    \b
        $ planemo clone --branch bwa-fix tools-iuc
        $ cd tools-iuc
        $ # Make changes.
        $ git add -p # Add desired changes.
        $ git commit -m "Fix bwa problem."
        $ planemo pull_request -m "Fix bwa problem."

    These changes do require that a github access token is
    specified in ~/.planemo.yml. An access token can be generated by going
    to https://github.com/settings/tokens.
    """
    github_util.pull_request(ctx, path, message=message, **kwds)
