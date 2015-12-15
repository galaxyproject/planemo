from ..deps import commands
from sys import platform as _platform
import os.path

IS_OS_X = _platform == "darwin"

# BSD 3-clause
CONDA_LICENSE = "http://docs.continuum.io/anaconda/eula"


def conda_link():
    if IS_OS_X:
        url = "https://repo.continuum.io/miniconda/Miniconda-latest-MacOSX-x86_64.sh"
    else:
        url = "https://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh"
    return url


def find_conda_prefix(conda_prefix=None):
    """ If supplied conda_prefix is not set, default to the default location
    for Miniconda installs.
    """
    if conda_prefix is None:
        return os.path.join(os.path.expanduser("~"), "miniconda2")
    return conda_prefix


def install_conda(prefix, shell_exec=None):
    if shell_exec is None:
        shell_exec = commands.shell
    download_cmd = " ".join(commands.download_command(conda_link(), quote=True))
    download_cmd = "%s > /tmp/conda.bash"
    install_cmd = "bash /tmp/conda.bash -b -p '%s'" % prefix
    full_command = "%s; %s" % (download_cmd, install_cmd)
    return shell_exec(full_command)
