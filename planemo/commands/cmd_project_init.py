"""Module describing the planemo ``project_init`` command."""
import os
import shutil
import tempfile

import click

from planemo import options
from planemo.cli import command_function
from planemo.io import (
    shell,
    untar_to,
    warn,
)

SOURCE_HOST = "https://codeload.github.com"
DOWNLOAD_URL = "%s/galaxyproject/planemo/tar.gz/master" % SOURCE_HOST
UNTAR_FILTER = "--strip-components=2"
UNTAR_ARGS = " -C %s -zxf - " + UNTAR_FILTER


@click.command("project_init")
@options.optional_project_arg(exists=None)
@click.option(
    '--template',
    default=None
)
@command_function
def cli(ctx, path, template=None, **kwds):
    """(Experimental) Initialize a new tool project.

    This is only a proof-of-concept demo right now.
    """
    if template is None:
        warn("Creating empty project, this function doesn't do much yet.")
    if not os.path.exists(path):
        os.makedirs(path)
    if template is None:
        return

    tempdir = tempfile.mkdtemp()
    try:
        untar_args = UNTAR_ARGS % (tempdir)
        untar_to(DOWNLOAD_URL, tempdir, untar_args)
        shell("ls '%s'" % (tempdir))
        shell("mv '%s/%s'/* '%s'" % (tempdir, template, path))
        shell("mv '%s/%s'/.* '%s'" % (tempdir, template, path))
    finally:
        shutil.rmtree(tempdir)
