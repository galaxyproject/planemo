"""
"""
import os
import tempfile
import shutil

import click

from planemo.cli import pass_context
from planemo import options
from planemo.io import (
    warn,
    untar_to,
    shell
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
@pass_context
def cli(ctx, path, template=None, **kwds):
    """Initialize a new tool project (demo only right now).
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
