"""
"""
import os
import shutil
import tempfile

import click

from planemo.cli import pass_context
from planemo.io import shell
from planemo import options
from planemo import shed


@click.command("shed_diff")
@options.optional_project_arg(exists=True)
@options.shed_owner_option()
@options.shed_name_option()
@options.shed_target_option()
@click.option(
    '--shed_target_source',
    help="Source Tool Shed to diff against (will ignore local project info"
         " specified). To compare the main Tool Shed against the test, set"
         " this to testtoolshed.",
    default=None,
)
@pass_context
def cli(ctx, path, **kwds):
    """Produce diff between local repository and Tool Shed contents.

    By default, this will produce a diff between this repository and what
    would be uploaded to the Tool Shed with the `shed_upload` command - but
    this command can be made to compare other combinations of repositories.
    Here are some examples::

        % # diff for this repository and the main Tool Shed
        % planemo shed_diff
        % # diff for this repository and the test Tool Shed
        % planemo shed_diff --shed_target testtoolshed
        % # diff for the test Tool Shed and main Tool Shed
        % planemo shed_diff --shed_target_source testtoolshed
        % # diff for two an explicitly specified repositories (ignores
        % # current project's shed YAML file.)
        % planemo shed_diff --owner peterjc --name blast_rbh
            --shed_target_source testtoolshed
    """
    working = tempfile.mkdtemp(prefix="tool_shed_diff_")
    try:
        diff_in(ctx, working, path, **kwds)
    finally:
        shutil.rmtree(working)


def diff_in(ctx, working, path, **kwds):
    shed_target_source = kwds.get("shed_target_source", None)

    label_a = "_%s_" % (shed_target_source if shed_target_source else "local")
    label_b = "_%s_" % kwds.get("shed_target", "B")

    mine = os.path.join(working, label_a)
    other = os.path.join(working, label_b)

    tsi = shed.tool_shed_client(ctx, read_only=True, **kwds)
    shed.download_tarball(
        ctx,
        tsi,
        path,
        destination=other,
        clean=True,
        **kwds
    )

    if shed_target_source:
        new_kwds = kwds.copy()
        new_kwds["shed_target"] = shed_target_source
        tsi = shed.tool_shed_client(ctx, read_only=True, **new_kwds)
        shed.download_tarball(
            ctx,
            tsi,
            path,
            destination=mine,
            clean=True,
            **new_kwds
        )
    else:
        tar_path = shed.build_tarball(path)
        cmd_template = 'mkdir "%s"; tar -xzf "%s" -C "%s"; rm -rf %s'
        shell(cmd_template % (mine, tar_path, mine, tar_path))

    shell('cd "%s"; diff -r %s %s' % (working, label_a, label_b))
