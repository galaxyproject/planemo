"""
"""
import os

import click

from planemo.cli import pass_context
from planemo import options
from planemo.io import warn, shell
from galaxy.tools.deps.commands import which

SOURCE_HOST = "https://codeload.github.com"
DOWNLOAD_URL = "%s/galaxyproject/planemo/tar.gz/master" % SOURCE_HOST
UNTAR_FILTER = "--strip-components=3 --wildcards --no-anchored '%s/**'"
UNTAR_CMD = "tar -C %s -zxvf - " + UNTAR_FILTER


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

    if which("wget"):
        download_cmd = "wget -O - %s"
    else:
        download_cmd = "curl %s"
    download_cmd = download_cmd % DOWNLOAD_URL
    untar_cmd = UNTAR_CMD % (path, template)
    shell("%s | %s" % (download_cmd, untar_cmd))
