"""Module describing the planemo ``shed_build`` command."""
import sys

import click
import shutil

from planemo.cli import command_function
from planemo import options
from planemo import shed


@click.command("shed_build")
@options.optional_tools_arg(multiple=False)
@command_function
def cli(ctx, path, **kwds):
    """Create a Galaxy tool tarball from a ``.shed.yml`` file.
    """

    def build(realized_repository):
        tarpath = shed.build_tarball(realized_repository.real_path)
        outpath = realized_repository.real_path + ".tar.gz"
        shutil.move(tarpath, outpath)
        print("Created: %s" % (outpath))
        return 0

    exit_code = shed.for_each_repository(ctx, build, [path], **kwds)
    sys.exit(exit_code)
