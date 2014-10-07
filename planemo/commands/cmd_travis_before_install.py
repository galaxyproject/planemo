import click

from planemo.cli import pass_context
from galaxy.tools.deps.commands import shell

FIX_EGGS_DIR = 'mkdir -p "$HOME/.python-eggs"; chmod 700 "$HOME/.python-eggs"'
CUSTOM_DEPS = (
    '[ -e ${TRAVIS_BUILD_DIR}/.travis/setup_custom_dependencies.bash ] && '
    '. ${TRAVIS_BUILD_DIR}/.travis/setup_custom_dependencies.bash'
)


@click.command('travis_before_install')
@pass_context
def cli(ctx):
    """This command is used internally by planemo to assist in contineous testing
    of tools with Travis CI (https://travis-ci.org/).
    """
    shell(FIX_EGGS_DIR)
    shell(CUSTOM_DEPS)
