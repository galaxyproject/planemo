"""Module describing the planemo ``conda_install`` command."""
import click

from galaxy.tools.deps import conda_util

from planemo import options
from planemo.cli import command_function
from planemo.conda import build_conda_context, collect_conda_targets
from planemo.exit_codes import EXIT_CODE_FAILED_DEPENDENCIES, ExitCodeException
from planemo.io import coalesce_return_codes, error


@click.command('conda_install')
@options.optional_tools_arg(multiple=True)
@options.conda_target_options()
@options.conda_auto_init_option()
@command_function
def cli(ctx, paths, **kwds):
    """Install conda packages for tool requirements."""
    conda_context = build_conda_context(ctx, **kwds)
    if not conda_context.is_conda_installed():
        auto_init = kwds.get("conda_auto_init", False)
        failed = True
        if auto_init:
            if conda_context.can_install_conda():
                if conda_util.install_conda(conda_context):
                    error("Attempted to install conda and failed.")
                else:
                    failed = False
            else:
                error("Cannot install conda, failing conda_install.")
        else:
            error("Conda not configured - run planemo conda_init' or pass --conda_auto_init to continue.")

        if failed:
            raise ExitCodeException(EXIT_CODE_FAILED_DEPENDENCIES)

    return_codes = []
    for conda_target in collect_conda_targets(ctx, paths):
        ctx.log("Install conda target %s" % conda_target)
        return_code = conda_util.install_conda_target(
            conda_target, conda_context=conda_context
        )
        return_codes.append(return_code)
    return coalesce_return_codes(return_codes, assert_at_least_one=True)
