"""Module describing the planemo ``dependency_script`` command."""
import os
import sys

from xml.etree import ElementTree as ET

import click

from planemo import options
from planemo.cli import command_function
from planemo.io import error, info
from planemo.shed2tap.base import BasePackage, Dependency


# We're using strict bash mode for dep_install.sh:
preamble_dep_install = """#!/bin/bash
#set -euo pipefail
set -eo pipefail
if [[ ! -d $INSTALL_DIR ]]
then
    echo "ERROR: Environment variable INSTALL_DIR not a directory!"
    exit 1
fi
# Set full strict mode now, side stepping case $INSTALL_DIR not setup.
set -euo pipefail
export DOWNLOAD_CACHE="${DOWNLOAD_CACHE:-./download_cache}"
if [[ ! -d $DOWNLOAD_CACHE ]]
then
    mkdir -p $DOWNLOAD_CACHE
fi
# Make this into an absolute path
export DOWNLOAD_CACHE=`(cd "$DOWNLOAD_CACHE"; pwd)`
echo "Using $DOWNLOAD_CACHE for cached downloads."
export INSTALL_DIR=`(cd "$INSTALL_DIR"; pwd)`
echo "Using $INSTALL_DIR for the installed files."
# Create a randomly named temp folder for working in
dep_install_tmp=${TMPDIR-/tmp}/dep_install.$RANDOM.$RANDOM.$RANDOM.$$
(umask 077 && mkdir $dep_install_tmp) || exit 1
"""

final_dep_install = """echo "Cleaning up..."
rm -rf $dep_install_tmp
echo "======================"
echo "Installation complete."
echo "======================"
"""

# Expect user to "source env.sh" so don't set strict mode,
# and don't use the exit command!
preamble_env_sh = """#!/bin/bash
if [[ ! -d $INSTALL_DIR ]]
then
    echo "ERROR: Environment variable INSTALL_DIR not a directory!"
fi
export INSTALL_DIR=${INSTALL_DIR:-$PWD}
export INSTALL_DIR=`(cd "$INSTALL_DIR"; pwd)`
"""


def find_tool_dependencis_xml(path, recursive):
    """Iterator function, quick & dirty tree walking."""
    if os.path.isfile(path):
        if os.path.basename(path) == "tool_dependencies.xml":
            yield path
    elif os.path.isdir(path):
        p = os.path.join(path, "tool_dependencies.xml")
        if os.path.isfile(p):
            yield p
        if recursive:
            for f in sorted(os.listdir(path)):
                p = os.path.join(path, f)
                if os.path.isdir(p):
                    # TODO: yield from
                    for x in find_tool_dependencis_xml(p, recursive):
                        yield x


def convert_tool_dep(dependencies_file):
    """Parse a tool_dependencies.xml into install.sh and env.sh commands.

    Returns two lists of strings, commands to add to install.sh and
    env.sh respectively.
    """
    install_cmds = []
    env_cmds = []

    root = ET.parse(dependencies_file).getroot()
    package_els = root.findall("package")

    packages = []
    dependencies = []
    for package_el in package_els:
        install_els = package_el.findall("install")
        assert len(install_els) in (0, 1)
        if len(install_els) == 0:
            repository_el = package_el.find("repository")
            assert repository_el is not None, "no repository in %s" % repository_el
            dependencies.append(Dependency(None, package_el, repository_el))
        else:
            install_el = install_els[0]
            packages.append(BasePackage(None, package_el, install_el, readme=None))

    if not packages:
        info("No packages in %s" % dependencies_file)
        return [], []

    assert len(packages) == 1, packages
    package = packages[0]
    name = package_el.attrib["name"]
    version = package_el.attrib["version"]

    # TODO - Set $INSTALL_DIR in the script
    # os.environ["INSTALL_DIR"] = os.path.abspath(os.curdir)
    for action in package.all_actions:
        inst, env = action.to_bash()
        install_cmds.extend(inst)
        env_cmds.extend(env)

    if install_cmds:
        install_cmds.insert(0, 'cd $dep_install_tmp')
        install_cmds.insert(0, 'specifc_action_done=0')
        install_cmds.insert(0, 'echo "%s"' % ('=' * 60))
        install_cmds.insert(0, 'echo "Installing %s version %s"' % (name, version))
        install_cmds.insert(0, 'echo "%s"' % ('=' * 60))
    if env_cmds:
        env_cmds.insert(0, 'specifc_action_done=0')
        env_cmds.insert(0, '#' + '=' * 60)
        env_cmds.insert(0, 'echo "Setting environment variables for %s version %s"' % (name, version))
        env_cmds.insert(0, '#' + '=' * 60)
        # TODO - define $INSTALL_DIR here?

    return install_cmds, env_cmds


def process_tool_dependencies_xml(tool_dep, install_handle, env_sh_handle):
    """Writes to handles, returns success as a boolean."""
    if not os.path.isfile(tool_dep):
        error('Missing file %s' % tool_dep)
        return False
    if not os.stat(tool_dep).st_size:
        error('Empty file %s' % tool_dep)
        return False
    try:
        install, env = convert_tool_dep(tool_dep)
    except Exception as err:
        # TODO - pass in ctx for logging?
        error('Error processing %s - %s' %
              (click.format_filename(tool_dep), err))
        if not isinstance(err, (NotImplementedError, RuntimeError)):
            # This is an unexpected error, traceback is useful
            import traceback
            error(traceback.format_exc() + "\n")
        return False
    # Worked...
    for cmd in install:
        install_handle.write(cmd + "\n")
    for cmd in env:
        env_sh_handle.write(cmd + "\n")
    return True


@click.command('dependency_script')
@options.shed_realization_options()
@options.dependencies_script_options()
@command_function
def cli(ctx, paths, recursive=False, fail_fast=True, download_cache=None):
    """Compile tool_dependencies.xml to bash script.

    An experimental approach parsing tool_dependencies.xml files into
    bash shell scripts, intended initially for use within Continuous
    Integration testing setups like TravisCI.

    Parses the ``tool_dependencies.xml`` files from the specified projects,
    and converts them into an installation bash script (``dep_install.sh``),
    and a shell script (``env.sh``) defining any new/edited environment
    variables.

    These are intended to be used via ``bash dep_install.sh`` (once), and as
    ``source env.sh`` prior to running any of the dependencies to set the
    environment variable within the current shell session.

    Both ``dep_install.sh`` and ``env.sh`` require ``$INSTALL_DIR`` be defined
    before running them, set to an existing directory with write permissions.
    Beware than if run on multiple tools, they can over-write each other (for
    example if you have packages for different versions of the same tool). In
    this case make separate calls to ``planemo dependency_script`` and call
    the scripts with different installation directories.

    This command will download (and cache) any URLs specified via Galaxy
    download actions. This is in order to decompress them and determine the
    relevant sub-folder to change into as per the Tool Shed install mechanism,
    so that this can be recorded as a ``cd`` comand in the bash script.

    The download cache used by ``planemo dependency_script`` and the resulting
    output script ``dep_install.sh`` defaults to ``./download_cache`` (under
    the current working directory), and can be set with ``$DOWNLOAD_CACHE``.

    If the ``tool_dependencies.xml`` file includes SHA256 checksums for
    downloads, these will be verified after downloading to the cache (by
    either ``planemo dependency_script`` or ``bash dep_install.sh``).

    This is experimental, and is initially intended for use within continuous
    integration testing setups like TravisCI to both verify the dependency
    installation receipe works, and to use this to run functional tests.
    """
    # TODO: Command line API for bash output filanames & install dir, cache.
    if download_cache:
        assert os.path.isdir(download_cache), download_cache
        # Effectively using this as a global variable, refactor this
        # once using a visitor pattern instead of action.to_bash()
        os.environ["DOWNLOAD_CACHE"] = os.path.abspath(download_cache)
        print("Using $DOWNLOAD_CACHE=%r" % os.environ["DOWNLOAD_CACHE"])
    failed = False
    with open("env.sh", "w") as env_sh_handle:
        with open("dep_install.sh", "w") as install_handle:
            install_handle.write(preamble_dep_install)
            env_sh_handle.write(preamble_env_sh)
            for path in paths:
                # ctx.log("Checking: %r" % path)
                if failed and fail_fast:
                    break
                for tool_dep in find_tool_dependencis_xml(path, recursive):
                    passed = process_tool_dependencies_xml(tool_dep,
                                                           install_handle,
                                                           env_sh_handle)
                    if passed:
                        info('Processed %s' % tool_dep)
                    else:
                        failed = True
                        if fail_fast:
                            for line in [
                                    '#' + '*' * 60,
                                    'echo "WARNING: Skipping %s"' % tool_dep,
                                    '#' + '*' * 60]:
                                install_handle.write(line + "\n")
                            break
                        # error("%s failed" % tool_dep)
            install_handle.write(final_dep_install)
    ctx.log("The End")
    if failed:
        error('Error processing one or more tool_dependencies.xml files.')
        sys.exit(1)
