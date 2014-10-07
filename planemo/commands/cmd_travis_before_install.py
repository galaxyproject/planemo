import os
import string

import click

from planemo.cli import pass_context
from planemo.io import shell

SETUP_FILE_NAME = "setup_custom_dependencies.bash"
SAMTOOLS_URL = (
    "http://archive.ubuntu.com/ubuntu/pool/universe/"
    "s/samtools/samtools_0.1.19-1_amd64.deb"
)

FIX_EGGS_DIR = 'mkdir -p "$HOME/.python-eggs"; chmod 700 "$HOME/.python-eggs"'
# samtools essentially required by Galaxy
INSTALL_SAMTOOLS = (
    "wget %s; "
    "sudo dpkg -i samtools_0.1.19-1_amd64.deb"
) % SAMTOOLS_URL

BUILD_ENVIRONMENT_TEMPLATE = """
export PATH=$PATH:${BUILD_BIN_DIR}
"""


@click.command('travis_before_install')
@pass_context
def cli(ctx):
    """This command is used internally by planemo to assist in contineous testing
    of tools with Travis CI (https://travis-ci.org/).
    """
    build_dir = os.environ.get("TRAVIS_BUILD_DIR", None)
    if not build_dir:
        raise Exception("Failed to determine ${TRAVIS_BUILD_DIR}")

    build_travis_dir = os.path.join(build_dir, ".travis")
    if not os.path.exists(build_travis_dir):
        os.makedirs(build_travis_dir)

    build_bin_dir = os.path.join(build_travis_dir, "bin")
    if not os.path.exists(build_bin_dir):
        os.makedirs(build_bin_dir)

    build_env_path = os.path.join(build_travis_dir, "env.sh")
    template_vars = {
        "BUILD_TRAVIS_DIR": build_travis_dir,
        "BUILD_BIN_DIR": build_bin_dir,
        "BUILD_ENV_PATH": build_env_path,
    }
    build_env = string.Template(BUILD_ENVIRONMENT_TEMPLATE).safe_substitute(
        **template_vars
    )
    open(build_env_path, "a").write(build_env)

    shell(FIX_EGGS_DIR)
    shell(INSTALL_SAMTOOLS)
    setup_file = os.path.join(build_travis_dir, SETUP_FILE_NAME)
    if os.path.exists(setup_file):
        shell(
            ". %s && bash -x %s" % (build_env_path, setup_file),
            env=template_vars
        )
