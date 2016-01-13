""" Click definitions for various shared options and arguments.
"""

import functools
import os

import click
from galaxy.tools.deps import docker_util
from .config import get_default_callback


def force_option(what="files"):
    return click.option(
        "-f",
        "--force",
        is_flag=True,
        help="Overwrite existing %s if present." % what,
    )


def skip_venv_option():
    return click.option(
        "--skip_venv",
        is_flag=True,
        help=("Do not create or source a virtualenv environment for Galaxy, "
              "this should be used or instance to preserve an externally "
              "configured virtual environment or conda environment.")
    )


def test_data_option():
    return click.option(
        "--test_data",
        type=click.Path(exists=True, file_okay=False, resolve_path=True),
        help="test-data directory to for specified tool(s).",
    )


def tool_data_table_option():
    return click.option(
        "--tool_data_table",
        type=click.Path(exists=True, file_okay=True, resolve_path=True),
        help="tool_data_table_conf.xml file to for specified tool(s).",
    )


def galaxy_root_option():
    return click.option(
        "--galaxy_root",
        type=click.Path(exists=True, file_okay=False, resolve_path=True),
        help="Root of development galaxy directory to execute command with.",
    )


def galaxy_sqlite_database_option():
    return click.option(
        "--galaxy_sqlite_database",
        callback=get_default_callback(None, resolve_path=True),
        type=click.Path(exists=True, file_okay=False, resolve_path=True),
        help="Preseeded Galaxy sqlite database to target.",
    )


def galaxy_cwl_root_option():
    return click.option(
        "--cwl_galaxy_root",
        type=click.Path(exists=True, file_okay=False, resolve_path=True),
        help=("Root of development galaxy directory to execute command with"
              " (must be branch of Galaxy with CWL support, this option"
              " is experimental and will be replaced with --galaxy_root when"
              " and if CWL support is merged into Galaxy."),
    )


def galaxy_port_option():
    return click.option(
        "--port",
        type=int,
        default=None,
        callback=get_default_callback("9090"),
        help="Port to serve Galaxy on (default is 9090).",
    )


def galaxy_host_option():
    return click.option(
        "--host",
        type=str,
        default=None,
        callback=get_default_callback("127.0.0.1"),
        help=("Host to bind Galaxy to. Default is 127.0.0.1 that is "
              "restricted to localhost connections for security reasons "
              "set to 0.0.0.0 to bind Galaxy to all ports including "
              "potentially publicly accessible ones."),
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


def enable_cwl_option():
    return click.option(
        "--cwl",
        is_flag=True,
        help=("Configure Galaxy for use with CWL tool."
              " (this option is experimental and will be replaced when"
              " and if CWL support is merged into Galaxy)."),
    )


def cwl_conformance_test():
    return click.option(
        "--conformance-test",
        is_flag=True,
        help=("Generate CWL conformance test object describing job. "
              "Required by CWL conformance test suite and implemented "
              "by cwltool reference implementation."),
    )


def brew_dependency_resolution():
    return click.option(
        "--brew_dependency_resolution",
        is_flag=True,
        help="Configure Galaxy to use plain brew dependency resolution.",
    )


def conda_dependency_resolution():
    return click.option(
        "--conda_dependency_resolution",
        is_flag=True,
        help="Configure Galaxy to use only conda for dependency resolution.",
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


def no_cache_galaxy_option():
    return click.option(
        "--no_cache_galaxy",
        is_flag=True,
        help=("Skip caching of Galaxy source and dependencies obtained with "
              "--install_galaxy. Not caching this results in faster "
              "downloads (no git) - so is better on throw away instances such "
              "with TravisCI. ")
    )


def galaxy_branch_option():
    return click.option(
        "--galaxy_branch",
        callback=get_default_callback(None),
        help=("Branch of Galaxy to target (defaults to master) if a Galaxy "
              "root isn't specified.")
    )


def galaxy_source_option():
    return click.option(
        "--galaxy_source",
        callback=get_default_callback(None),
        help=("Git source of Galaxy to target (defaults to the official "
              "galaxyproject github source if a Galaxy root isn't "
              "specified.")
    )


def skip_install_option():
    return click.option(
        "--skip_install",
        is_flag=True,
        help="Skip installation - only source requirements already available."
    )


def brew_option():
    return click.option(
        "--brew",
        type=click.Path(exists=True, file_okay=True, dir_okay=False),
        help="Homebrew 'brew' executable to use."
    )


def conda_prefix_option():
    return click.option(
        "--conda_prefix",
        type=click.Path(file_okay=False, dir_okay=True),
        help="Conda prefix to use for conda dependency commands."
    )


def conda_exec_option():
    return click.option(
        "--conda_exec",
        type=click.Path(exists=True, file_okay=True, dir_okay=False),
        help="Location of conda executable."
    )


def conda_debug_option():
    return click.option(
        "--conda_debug",
        is_flag=True,
        help="Enable more verbose conda logging."
    )


def conda_ensure_channels_option():
    return click.option(
        "--conda_ensure_channels",
        type=str,
        help=("Ensure conda is configured with specified comma separated "
              "list of channels."),
        default="r,bioconda"
    )


def conda_copy_dependencies_option():
    return click.option(
        "--conda_copy_dependencies",
        is_flag=True,
        help=("Conda dependency resolution for Galaxy will copy dependencies "
              "instead of attempting to link them.")
    )


def conda_auto_install_option():
    return click.option(
        "--conda_auto_install",
        is_flag=True,
        help=("Conda dependency resolution for Galaxy will auto install "
              "will attempt to install requested but missing packages.")
    )


def conda_auto_init_option():
    return click.option(
        "--conda_auto_init",
        is_flag=True,
        help=("Conda dependency resolution for Galaxy will auto install "
              "conda itself using miniconda if not availabe on conda_prefix.")
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
    return click.argument("path", metavar="TOOL_PATH", type=arg_type)


def required_job_arg():
    """ Decorate click method as requiring the path to a single tool.
    """
    arg_type = click.Path(
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
    )
    return click.argument("job_path", metavar="JOB_PATH", type=arg_type)


def _optional_tools_default(ctx, param, value):
    if param.name == "paths" and len(value) == 0:
        return [os.path.abspath(os.getcwd())]
    else:
        return value


def optional_tools_arg(multiple=False):
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
    name = "paths" if multiple else "path"
    nargs = -1 if multiple else 1
    return click.argument(
        name,
        metavar="TOOL_PATH",
        type=arg_type,
        nargs=nargs,
        callback=_optional_tools_default,
    )


class ProjectOrRepositry(click.Path):

    def __init__(self, **kwds):
        super(ProjectOrRepositry, self).__init__(**kwds)

    def convert(self, value, param, ctx):
        if value and value.startswith("git:") or value.startswith("git+"):
            return value
        else:
            return super(ProjectOrRepositry, self).convert(value, param, ctx)


def shed_project_arg(multiple=True):
    arg_type = ProjectOrRepositry(
        exists=True,
        file_okay=False,
        dir_okay=True,
        writable=True,
        resolve_path=True,
    )
    name = "paths" if multiple else "path"
    nargs = -1 if multiple else 1
    return click.argument(
        name,
        metavar="PROJECT",
        type=arg_type,
        nargs=nargs,
        callback=_optional_tools_default,
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
        "path",
        metavar="PROJECT",
        default=".",
        type=arg_type
    )


def no_cleanup_option():
    return click.option(
        "--no_cleanup",
        is_flag=True,
        help=("Do not cleanup temp files created for and by Galaxy.")
    )


def docker_cmd_option():
    return click.option(
        "--docker_cmd",
        default=docker_util.DEFAULT_DOCKER_COMMAND,
        help="Command used to launch docker (defaults to docker)."
    )


def docker_sudo_option():
    return click.option(
        "--docker_sudo",
        is_flag=True,
        help="Flag to use sudo when running docker."
    )


def docker_sudo_cmd_option():
    return click.option(
        "--docker_sudo_cmd",
        default=None,
        help="sudo command to use when --docker_sudo is enabled " +
             "(defaults to sudo).",
        callback=get_default_callback(docker_util.DEFAULT_SUDO_COMMAND),
    )


def docker_host_option():
    return click.option(
        "--docker_host",
        default=None,
        help="Docker host to target when executing docker commands " +
             "(defaults to localhost).",
        callback=get_default_callback(docker_util.DEFAULT_HOST),
    )


def shed_owner_option():
    return click.option(
        "--owner",
        help="Tool Shed repository owner (username)."
    )


def shed_name_option():
    return click.option(
        "--name",
        help="Tool Shed repository name (defaults to the inferred "
             "tool directory name)."
    )


def validate_shed_target_callback(ctx, param, value):
    target = get_default_callback("toolshed")(ctx, param, value)
    if target is None:
        ctx.fail("default_shed_target set to None, must specify a value for "
                 "--shed_target to run this command.")
    return target


def shed_target_option():
    return click.option(
        "-t",
        "--shed_target",
        help="Tool Shed to target (this can be 'toolshed', 'testtoolshed', "
             "'local' (alias for http://localhost:9009/), an arbitrary url "
             "or mappings defined ~/.planemo.yml.",
        default=None,
        callback=validate_shed_target_callback,
    )


def shed_key_option():
    return click.option(
        "--shed_key",
        help=("API key for Tool Shed access. An API key is required unless "
              "e-mail and password is specified. This key can be specified "
              "with either --shed_key or --shed_key_from_env.")
    )


def shed_key_from_env_option():
    return click.option(
        "--shed_key_from_env",
        help="Environment variable to read API key for Tool Shed access from."
    )


def shed_email_option():
    return click.option(
        "--shed_email",
        help="E-mail for Tool Shed auth (required unless shed_key is "
             "specified)."
    )


def shed_password_option():
    return click.option(
        "--shed_password",
        help="Password for Tool Shed auth (required unless shed_key is "
             "specified)."
    )


def shed_skip_upload():
    return click.option(
        "--skip_upload",
        is_flag=True,
        help=("Skip upload contents as part of operation, only update "
              "metadata.")
    )


def shed_skip_metadata():
    return click.option(
        "--skip_metadata",
        is_flag=True,
        help=("Skip metadata update as part of operation, only upload "
              "new contents.")
    )


def shed_message_option():
    return click.option(
        "-m",
        "--message",
        help="Commit message for tool shed upload."
    )


def shed_force_create_option():
    return click.option(
        "--force_repository_creation",
        help=("If a repository cannot be found for the specified user/repo "
              "name pair, then automatically create the repository in the "
              "toolshed."),
        is_flag=True,
        default=False
    )


def shed_check_diff_option():
    return click.option(
        "--check_diff",
        is_flag=True,
        help=("Skip uploading if the shed_diff detects there would be no "
              "'difference' (only attributes populated by the shed would "
              "be updated.)")
    )


def shed_upload_options():
    return _compose(
        shed_message_option(),
        shed_force_create_option(),
        shed_check_diff_option(),
    )


def shed_realization_options():
    return _compose(
        shed_project_arg(multiple=True),
        recursive_shed_option(),
        shed_fail_fast_option(),
    )


def shed_repo_options():
    return _compose(
        shed_owner_option(),
        shed_name_option(),
    )


def shed_publish_options():
    """ Common options for commands that require publishing to a
    a shed.
    """
    return _compose(
        shed_realization_options(),
        shed_repo_options(),
        shed_target_options(),
    )


def shed_read_options():
    """ Common options that require read access to mapped repositories
    in a shed.
    """
    return _compose(
        shed_realization_options(),
        shed_repo_options(),
        shed_target_options(),
    )


def shed_target_options():
    """ Common options for commands that require read-only
    interactions with a shed.
    """
    return _compose(
        shed_email_option(),
        shed_key_option(),
        shed_key_from_env_option(),
        shed_password_option(),
        shed_target_option(),
    )


def conda_target_options():
    return _compose(
        conda_prefix_option(),
        conda_exec_option(),
        conda_debug_option(),
        conda_ensure_channels_option(),
    )


def galaxy_run_options():
    return _compose(
        galaxy_target_options(),
        galaxy_port_option(),
        galaxy_host_option(),
    )


def galaxy_config_options():
    return _compose(
        test_data_option(),
        tool_data_table_option(),
        dependency_resolvers_option(),
        tool_dependency_dir_option(),
        brew_dependency_resolution(),
        shed_dependency_resolution(),
        conda_target_options(),
        conda_dependency_resolution(),
        conda_copy_dependencies_option(),
        conda_auto_install_option(),
        conda_auto_init_option(),
    )


def galaxy_target_options():
    return _compose(
        galaxy_root_option(),
        galaxy_sqlite_database_option(),
        install_galaxy_option(),
        galaxy_branch_option(),
        galaxy_source_option(),
        skip_venv_option(),
        no_cache_galaxy_option(),
        no_cleanup_option(),
        job_config_option(),
    )


def galaxy_serve_options():
    return _compose(
        galaxy_run_options(),
        galaxy_config_options(),
    )


def shed_fail_fast_option():
    return click.option(
        "--fail_fast",
        is_flag=True,
        default=False,
        help="If multiple repositories are specified and an error occurs "
             "stop immediately instead of processing remaining repositories."
    )


def lint_xsd_option():
    return click.option(
        "--xsd",
        is_flag=True,
        default=False,
        help=("Include experimental tool XSD validation in linting "
              "process (requires xmllint on PATH or lxml installed).")
    )


def report_level_option():
    return click.option(
        "--report_level",
        type=click.Choice(["all", "warn", "error"]),
        default="all",
    )


def report_xunit():
    return click.option(
        "--report_xunit",
        type=click.Path(file_okay=True, resolve_path=True),
        help="Output an XUnit report, useful for CI testing",
        default=None,
    )


def skip_option():
    return click.option(
        "-s",
        "--skip",
        default=None,
        help=("Comma-separated list of lint tests to skip (e.g send ."
              "--skip 'citations,xml_order' to skip linting of citations "
              "and best-practice XML ordering.")
    )


def fail_level_option():
    return click.option(
        "--fail_level",
        type=click.Choice(['warn', 'error']),
        default="warn"
    )


def recursive_shed_option():
    return recursive_option(
        "Recursively perform command for nested repository directories.",
    )


def recursive_option(help="Recursively perform command for subdirectories."):
    return click.option(
        "-r",
        "--recursive",
        is_flag=True,
        help=help,
    )


def tool_test_json():
    target_path = click.Path(
        file_okay=True,
        dir_okay=False,
        resolve_path=True,
    )
    return click.argument(
        'path',
        metavar="FILE_PATH",
        type=target_path,
        default="tool_test_output.json",
    )


def test_report_options():
    return _compose(
        click.option(
            "--test_output",
            type=click.Path(file_okay=True, resolve_path=True),
            callback=get_default_callback("tool_test_output.html",
                                          resolve_path=True),
            help=("Output test report (HTML - for humans) defaults to "
                  "tool_test_output.html."),
            default=None,
        ),
        click.option(
            "--test_output_text",
            type=click.Path(file_okay=True, resolve_path=True),
            callback=get_default_callback(None, resolve_path=True),
            help=("Output test report (Basic text - for display in CI)"),
            default=None,
        ),
        click.option(
            "--test_output_markdown",
            type=click.Path(file_okay=True, resolve_path=True),
            callback=get_default_callback(None, resolve_path=True),
            help=("Output test report (Markdown style - for humans & "
                  "computers)"),
            default=None,
        ),
    )


def test_options():
    return _compose(
        click.option(
            "--update_test_data",
            is_flag=True,
            help="Update test-data directory with job outputs (normally"
                 " written to directory --job_output_files if specified.)"
        ),
        test_report_options(),
        click.option(
            "--test_output_xunit",
            type=click.Path(file_okay=True, resolve_path=True),
            callback=get_default_callback(None, resolve_path=True),
            help="Output test report (xUnit style - for computers).",
            default=None,
        ),
        click.option(
            "--test_output_json",
            type=click.Path(file_okay=True, resolve_path=True),
            callback=get_default_callback("tool_test_output.json",
                                          resolve_path=True),
            help=("Output test report (planemo json) defaults to "
                  "tool_test_output.json."),
            default=None,
        ),
        click.option(
            "--job_output_files",
            type=click.Path(file_okay=False, resolve_path=True),
            help="Write job outputs to specified directory.",
            default=None,
        ),
        click.option(
            "--summary",
            type=click.Choice(["none", "minimal", "compact"]),
            default="minimal",
            help=("Summary style printed to planemo's standard output (see "
                  "output reports for more complete summary). Set to 'none' "
                  "to disable completely.")
        )
    )


def _compose(*functions):
    def compose2(f, g):
        return lambda x: f(g(x))
    return functools.reduce(compose2, functions)


def dependencies_script_options():
    return _compose(
        click.option(
            "--download_cache",
            type=click.Path(file_okay=False, resolve_path=True),
            callback=get_default_callback(None),
            help=("Directory to cache downloaded files, default is $DOWNLOAD_CACHE"),
            default=None,
        ),
    )
