"""
"""
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
@options.optional_project_arg(exists=True)
@click.option(
    '--destination',
    default="shed_download.tar.gz",
    type=target_path,
    help="Destination of tarball to download - if this doesn't end in 'gz' it "
         "will be treated as a directory to extract tool contents into"
         "(defaults to shed_download.tar.gz)."
)
@options.shed_owner_option()
@options.shed_name_option()
@options.shed_target_option()
@pass_context
def cli(ctx, path, **kwds):
    """Download a tool repository as a tarball from the tool shed and extract
    to the specified directory.
    """
    tsi = shed.tool_shed_client(ctx, read_only=True, **kwds)
    shed.download_tarball(ctx, tsi, path, **kwds)
