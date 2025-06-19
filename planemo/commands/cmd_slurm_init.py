"""Module describing the planemo ``slurm_init`` command."""

import shutil
import tempfile

import click

from planemo.cli import command_function
from planemo.io import (
    info,
    shell,
    untar_to,
)

SLRUM_DRMAA_VERSION = "1.1.4"
DOWNLOAD_URL = f"https://github.com/natefoo/slurm-drmaa/releases/download/{SLRUM_DRMAA_VERSION}/slurm-drmaa-{SLRUM_DRMAA_VERSION}.tar.gz"


@click.command("slurm_init")
@command_function
def cli(ctx, **kwds):
    """Initialize a copy of the SLURM DRMAA library."""
    dest = f"{ctx.workspace}/libdrmaa.so"
    tempdir = tempfile.mkdtemp()
    tar_args = ["-zxf", "-", "--strip-components=1"]
    try:
        # This used to show the wget command (and in fact use wget), it doesn't anymore but
        # this should be restored. The point is to show developers what is happening and how to do it.
        untar_to(DOWNLOAD_URL, tar_args=tar_args, dest_dir=tempdir)
        shell(["mkdir", "dest"], cwd=tempdir)
        shell(["./autogen.sh"], cwd=tempdir)
        shell(["./configure", f"--prefix={tempdir}/dist"], cwd=tempdir)
        shell(["make"], cwd=tempdir)
        shell(["make install"], cwd=tempdir)
        shutil.move(f"{tempdir}/dist/lib/libdrmaa.so", dest)
        info(f"SLURM DRMAA library initialized successfully and copied to {dest}.")
    finally:
        shutil.rmtree(tempdir)
