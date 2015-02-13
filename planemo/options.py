""" Click definitions for various shared options and arguments.
"""
import click
from galaxy.tools.deps import docker_util


def test_data_option():
    return click.option(
        "--test_data",
        type=click.Path(exists=True, file_okay=False, resolve_path=True),
        help='test-data directory to for specified tool(s).'
    )


def tool_data_table_option():
    return click.option(
        "--tool_data_table",
        type=click.Path(exists=True, file_okay=True, resolve_path=True),
        help='tool_data_table_conf.xml file to for specified tool(s).'
    )


def galaxy_root_option():
    return click.option(
        "--galaxy_root",
        type=click.Path(exists=True, file_okay=False, resolve_path=True),
        help='Root of development galaxy directory to execute command with.'
    )


def dependency_resolvers_option():
    return click.option(
        "--dependency_resolvers_config_file",
        type=click.Path(
            exists=True,
            file_okay=True,
            dir_okay=False,
            resolve_path=True
        ),
        help="Dependency resolver configuration for Galaxy to target.",
    )


def brew_dependency_resolution():
    return click.option(
        "--brew_dependency_resolution",
        is_flag=True,
        help="Configure Galaxy to use plain brew dependency resolution.",
    )


def shed_dependency_resolution():
    return click.option(
        "--shed_dependency_resolution",
        is_flag=True,
        help=("Configure Galaxy to use brewed Tool Shed dependency"
              " resolution."),
    )


def job_config_option():
    return click.option(
        "--job_config_file",
        type=click.Path(
            exists=True,
            file_okay=True,
            dir_okay=False,
            resolve_path=True
        ),
        help="Job configuration file for Galaxy to target.",
    )


def tool_dependency_dir_option():
    return click.option(
        "--tool_dependency_dir",
        type=click.Path(
            exists=True,
            file_okay=False,
            dir_okay=True,
            resolve_path=True
        ),
        help="Tool dependency dir for Galaxy to target.",
    )


def install_galaxy_option():
    return click.option(
        "--install_galaxy",
        is_flag=True,
        help="Download and configure a disposable copy of Galaxy from github."
    )


def brew_option():
    return click.option(
        '--brew',
        type=click.Path(exists=True, file_okay=True, dir_okay=False),
        help="Homebrew 'brew' executable to use."
    )


def required_tool_arg():
    """ Decorate click method as requiring the path to a single tool.
    """
    arg_type = click.Path(
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
    )
    return click.argument('path', metavar="TOOL_PATH", type=arg_type)


def optional_tools_arg():
    """ Decorate click method as optionally taking in the path to a tool
    or directory of tools. If no such argument is given the current working
    directory will be treated as a directory of tools.
    """
    arg_type = click.Path(
        exists=True,
        file_okay=True,
        dir_okay=True,
        readable=True,
        resolve_path=True,
    )
    return click.argument(
        'path',
        metavar="TOOL_PATH",
        default=".",
        type=arg_type
    )


def optional_project_arg(exists=True):
    arg_type = click.Path(
        exists=exists,
        file_okay=False,
        dir_okay=True,
        writable=True,
        resolve_path=True,
    )
    return click.argument(
        'path',
        metavar="PROJECT",
        default=".",
        type=arg_type
    )


def docker_cmd_option():
    return click.option(
        '--docker_cmd',
        default=docker_util.DEFAULT_DOCKER_COMMAND,
        help="Command used to launch docker (defaults to docker)."
    )


def docker_sudo_option():
    return click.option(
        '--docker_sudo',
        is_flag=True,
        help="Flag to use sudo when running docker."
    )


def docker_sudo_cmd_option():
    return click.option(
        '--docker_sudo_cmd',
        default=docker_util.DEFAULT_SUDO_COMMAND,
        help="sudo command to use when --docker_sudo is enabled " +
             "(defaults to sudo)."
    )


def docker_host_option():
    return click.option(
        '--docker_host',
        default=docker_util.DEFAULT_HOST,
        help="Docker host to target when executing docker commands " +
             "(defaults to localhost)."
    )


def shed_owner_option():
    return click.option(
        '--owner',
        help="Tool Shed repository owner (username)."
    )


def shed_name_option():
    return click.option(
        '--name',
        help="Tool Shed repository name (defaults to the inferred "
             "tool directory name)."
    )


def shed_target_option():
    return click.option(
        '--shed_target',
        help="Tool Shed to target (this can be 'toolshed', 'testtoolshed', "
             "'local' (alias for http://localhost:9009/) or an arbitrary"
             "url).",
        default="toolshed",
    )
