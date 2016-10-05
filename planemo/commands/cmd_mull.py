"""Module describing the planemo ``mull`` command."""
import click

from galaxy.tools.deps.mulled.mulled_build import mull_targets
from galaxy.tools.deps.mulled.util import build_target

from planemo import options
from planemo.cli import command_function
from planemo.conda import collect_conda_target_lists
from planemo.mulled import build_mull_target_kwds


@click.command('mull')
@options.optional_tools_arg(multiple=True)
@options.recursive_option()
@command_function
def cli(ctx, paths, **kwds):
    """Build containers for specified tools.

    Supplied tools will be inspected for referenced requirement packages. For
    each combination of requirements a "mulled" container will be built. Galaxy
    can automatically discover this container and subsequently use it to run
    or test the tool.

    For this to work, the tool's requirements will need to be present in a known
    Conda channel such as bioconda (https://github.com/bioconda/bioconda-recipes).
    This can be verified by running ``planemo lint --conda_requirements`` on the
    target tool(s).
    """
    for conda_targets in collect_conda_target_lists(ctx, paths):
        mulled_targets = map(lambda c: build_target(c.package, c.version), conda_targets)
        mull_target_kwds = build_mull_target_kwds(ctx, **kwds)
        mull_targets(mulled_targets, command="build", **mull_target_kwds)
