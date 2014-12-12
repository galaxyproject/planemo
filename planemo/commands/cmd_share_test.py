"""
"""
import click

from planemo.cli import pass_context
from planemo.io import info
from planemo import github_util

PLANEMO_TEST_VIEWER_URL_TEMPLATE = (
    "http://galaxyproject.github.io/planemo/tool_test_viewer.html"
    "?test_data_url=%s"
)

target_path = click.Path(
    file_okay=True,
    dir_okay=False,
    resolve_path=True,
)


@click.command("share_test")
@click.argument(
    'path',
    metavar="FILE_PATH",
    type=target_path,
    default="tool_test_output.json",
)
@pass_context
def cli(ctx, path, **kwds):
    """Publish JSON test results to Github Gist and produce sharable URL.

    Sharable URL  can be used to share an HTML version of the report that
    can be easily embedded in pull requests or commit messages. Requires
    a ~/.planemo.yml with Github 'username' and 'password' defined in a
    'github' section of that configuration file.
    """
    file_url = github_util.publish_as_gist_file(ctx, path)
    share_url = PLANEMO_TEST_VIEWER_URL_TEMPLATE % file_url
    info("File published to Github Gist.")
    info("Raw URL: %s" % file_url)
    info("Share results with URL: %s" % share_url)
    markdown = "[View Tool Test Results](%s)" % share_url
    info("Embed results with markdown: %s" % markdown)
