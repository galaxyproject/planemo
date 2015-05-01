"""
"""
import sys

import click

from planemo.cli import pass_context
from planemo import options
from planemo import shed

target_path = click.Path(
    file_okay=True,
    writable=True,
    resolve_path=True,
)


@click.command("shed_download")
@options.shed_project_arg()
@click.option(
    '--destination',
    default="shed_download.tar.gz",
    type=target_path,
    help="Destination pattern of tarball(s) to download - if this doesn't "
         "end in 'gz' it will be treated as a directory to extract tool "
         "contents into (defaults to shed_download.tar.gz). If multiple "
         "repositories are discovered in a .shed.yml file these will be "
         "created as shed_download_<name>.tar.gz by default for instance, "
         "simpler repositories will just be downloaded to the specified file."
)
@options.shed_owner_option()
@options.shed_name_option()
@options.shed_target_option()
@options.recursive_shed_option()
@options.shed_fail_fast_option()
@pass_context
def cli(ctx, path, **kwds):
    """Download a tool repository as a tarball from the tool shed and extract
    to the specified directory.
    """
    tsi = shed.tool_shed_client(ctx, read_only=True, **kwds)

    def download(realized_repository):
        return shed.download_tarball(ctx, tsi, realized_repository, **kwds)

    exit_code = shed.for_each_repository(download, path, **kwds)
    sys.exit(exit_code)
