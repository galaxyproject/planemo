"""Module describing the planemo ``conda_install`` command."""
import click

from planemo.cli import command_function
from planemo.io import coalesce_return_codes
from planemo import options

from planemo.conda import build_conda_context, collect_conda_targets

from galaxy.tools.deps import conda_util


@click.command('conda_install')
@options.optional_tools_arg()
@options.conda_target_options()
@command_function
def cli(ctx, path, **kwds):
    """Install conda packages for tool requirements."""
    conda_context = build_conda_context(**kwds)
    return_codes = []
    for conda_target in collect_conda_targets(path):
        ctx.log("Install conda target %s" % conda_target)
        return_code = conda_util.install_conda_target(
            conda_target, conda_context=conda_context
        )
        return_codes.append(return_code)
    return coalesce_return_codes(return_codes, assert_at_least_one=True)
