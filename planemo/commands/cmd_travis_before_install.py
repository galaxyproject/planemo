"""Module describing the planemo ``travis_before_install`` command."""
import os
import string

import click

from planemo.cli import command_function
from planemo.io import shell

SETUP_FILE_NAME = "setup_custom_dependencies.bash"
SAMTOOLS_DEB = 'samtools_0.1.19-1_amd64.deb'
SAMTOOLS_URL = "http://archive.ubuntu.com/ubuntu/pool/universe/s/samtools/%s" % SAMTOOLS_DEB

BUILD_ENVIRONMENT_TEMPLATE = """
export PATH=$PATH:${BUILD_BIN_DIR}
"""


@click.command('travis_before_install')
@command_function
def cli(ctx):
    """Internal command for GitHub/TravisCI testing.

    This command is used internally by planemo to assist in continuous testing
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
    with open(build_env_path, "a") as fh:
        fh.write(build_env)

    eggs_dir = os.path.join(os.getenv('HOME'), '.python-eggs')
    if not os.path.exists(eggs_dir):
        os.makedirs(eggs_dir, 0o700)
    else:
        os.chmod(eggs_dir, 0o700)
    # samtools essentially required by Galaxy
    shell(['wget', SAMTOOLS_URL])
    shell(['sudo', 'dpkg', '-i', SAMTOOLS_DEB])
    setup_file = os.path.join(build_travis_dir, SETUP_FILE_NAME)
    if os.path.exists(setup_file):
        env = template_vars
        env['PATH'] = build_bin_dir
        shell(['bash', '-x', setup_file], env=env)
