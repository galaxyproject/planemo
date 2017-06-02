"""Click definitions for various shared options and arguments."""

from __future__ import absolute_import

import functools
import os

import click

from galaxy.tools.deps import docker_util

from .config import planemo_option


def force_option(what="files"):
    return planemo_option(
        "-f",
        "--force",
        is_flag=True,
        help="Overwrite existing %s if present." % what,
    )


def skip_venv_option():
    return planemo_option(
        "--skip_venv",
        is_flag=True,
        help=("Do not create or source a virtualenv environment for Galaxy, "
              "this should be used or instance to preserve an externally "
              "configured virtual environment or conda environment.")
    )


def run_engine_option():
    return planemo_option(
        "--engine",
        type=click.Choice(["galaxy", "docker_galaxy", "cwltool"]),
        default="galaxy",
        use_global_config=True,
        help=("Select an engine to run or test aritfacts such as tools "
              "and workflows. Defaults to a local Galaxy, but running Galaxy within "
              "a Docker container or the CWL reference implementation 'cwltool' and "
              "be selected.")
    )


def non_strict_cwl_option():
    return planemo_option(
        "--non_strict_cwl",
        default=False,
        is_flag=True,
        help="Disable strict validation of CWL.",
    )


def serve_engine_option():
    return planemo_option(
        "--engine",
        type=click.Choice(["galaxy", "docker_galaxy"]),
        default="galaxy",
        use_global_config=True,
        use_env_var=True,
        help=("Select an engine to serve aritfacts such as tools "
              "and workflows. Defaults to a local Galaxy, but running Galaxy within "
              "a Docker container.")
    )


def cwltool_no_container_option():
    return planemo_option(
        "--no-container",
        "--no_container",
        is_flag=True,
        default=False,
        use_global_config=True,
        help=("If cwltool engine is used, disable Docker container usage.")
    )


def test_data_option():
    return planemo_option(
        "--test_data",
        type=click.Path(exists=True, file_okay=False, resolve_path=True),
        help="test-data directory to for specified tool(s).",
    )


def extra_tools_option():
    return planemo_option(
        "--extra_tools",
        type=click.Path(exists=True,
                        file_okay=True,
                        dir_okay=True,
                        resolve_path=True),
        multiple=True,
        help=("Extra tool sources to include in Galaxy's tool panel (file or "
              "directory). These will not be linted/tested/etc... but they "
              "will be available to workflows and for interactive use.")
    )


def tool_data_table_option():
    return planemo_option(
        "--tool_data_table",
        type=click.Path(exists=True, file_okay=True, resolve_path=True),
        help="tool_data_table_conf.xml file to for specified tool(s).",
    )


def galaxy_email_option():
    return planemo_option(
        "--galaxy_email",
        type=str,
        default="planemo@galaxyproject.org",
        use_global_config=True,
        use_env_var=True,
        help="E-mail address to use when launching single-user Galaxy server.",
    )


def galaxy_root_option():
    return planemo_option(
        "--galaxy_root",
        type=click.Path(exists=True, file_okay=False, resolve_path=True),
        help="Root of development galaxy directory to execute command with.",
    )


def galaxy_database_seed_option():
    return planemo_option(
        "--galaxy_database_seed",
        default=None,
        use_global_config=True,
        use_env_var=True,
        type=click.Path(exists=True, file_okay=True, resolve_path=True),
        help="Preseeded Galaxy sqlite database to target.",
    )


def galaxy_cwl_root_option():
    return planemo_option(
        "--cwl_galaxy_root",
        use_global_config=True,
        type=click.Path(exists=True, file_okay=False, resolve_path=True),
        help=("Root of development galaxy directory to execute command with"
              " (must be branch of Galaxy with CWL support, this option"
              " is experimental and will be replaced with --galaxy_root when"
              " and if CWL support is merged into Galaxy."),
    )


def galaxy_port_option():
    return planemo_option(
        "--port",
        type=int,
        default="9090",
        use_global_config=True,
        help="Port to serve Galaxy on (default is 9090).",
    )


def galaxy_host_option():
    return planemo_option(
        "--host",
        type=str,
        default="127.0.0.1",
        use_global_config=True,
        help=("Host to bind Galaxy to. Default is 127.0.0.1 that is "
              "restricted to localhost connections for security reasons "
              "set to 0.0.0.0 to bind Galaxy to all ports including "
              "potentially publicly accessible ones."),
    )


def dependency_resolvers_option():
    return planemo_option(
        "--dependency_resolvers_config_file",
        type=click.Path(
            exists=True,
            file_okay=True,
            dir_okay=False,
            resolve_path=True
        ),
        use_global_config=True,
        help="Dependency resolver configuration for Galaxy to target.",
    )


def enable_cwl_option():
    return planemo_option(
        "--cwl",
        is_flag=True,
        help=("Configure Galaxy for use with CWL tool."
              " (this option is experimental and will be replaced when"
              " and if CWL support is merged into Galaxy)."),
    )


def build_cwl_option():
    return planemo_option(
        "--cwl",
        is_flag=True,
        help="Build a CWL tool instead of a Galaxy tool.",
    )


def run_output_directory_option():
    return planemo_option(
        "output_directory",
        "--output_directory",
        "--outdir",
        type=click.Path(
            file_okay=False,
            dir_okay=True,
            resolve_path=True,
        ),
        default=None,
        help=("Where to store outputs of a 'run' task."),
    )


def run_output_json_option():
    return planemo_option(
        "output_json",
        "--output_json",
        type=click.Path(
            file_okay=True,
            dir_okay=False,
            resolve_path=True,
        ),
        default=None,
        help=("Where to store JSON dictionary describing outputs of "
              "a 'run' task."),
    )


def no_dependency_resolution():
    return planemo_option(
        "--no_dependency_resolution",
        is_flag=True,
        help="Configure Galaxy with no dependency resolvers.",
    )


def brew_dependency_resolution():
    return planemo_option(
        "--brew_dependency_resolution",
        is_flag=True,
        help="Configure Galaxy to use plain brew dependency resolution.",
    )


def conda_dependency_resolution():
    return planemo_option(
        "--conda_dependency_resolution",
        is_flag=True,
        help="Configure Galaxy to use only conda for dependency resolution.",
    )


def shed_dependency_resolution():
    return planemo_option(
        "--shed_dependency_resolution",
        is_flag=True,
        help=("Configure Galaxy to use brewed Tool Shed dependency"
              " resolution."),
    )


def file_path_option():
    return planemo_option(
        "--file_path",
        type=click.Path(
            file_okay=False,
            dir_okay=True,
            resolve_path=True
        ),
        help="Location for files created by Galaxy (e.g. database/files).",
        default=None,
        use_global_config=True,
    )


def database_connection_option():
    return planemo_option(
        "--database_connection",
        type=str,
        help="Database connection string to use for Galaxy.",
        default=None,
        use_global_config=True,
    )


def shed_tools_conf_option():
    return planemo_option(
        "--shed_tool_conf",
        type=str,
        help="Location of shed tools conf file for Galaxy.",
        default=None,
        use_global_config=True,
    )


def shed_tools_directory_option():
    return planemo_option(
        "--shed_tool_path",
        type=str,
        help="Location of shed tools directory for Galaxy.",
        default=None,
        use_global_config=True,
    )


def tool_dependency_dir_option():
    return planemo_option(
        "--tool_dependency_dir",
        type=click.Path(
            exists=True,
            file_okay=False,
            dir_okay=True,
            resolve_path=True
        ),
        default=None,
        use_global_config=True,
        help="Tool dependency dir for Galaxy to target.",
    )


def job_config_option():
    return planemo_option(
        "--job_config_file",
        type=click.Path(
            exists=True,
            file_okay=True,
            dir_okay=False,
            resolve_path=True
        ),
        help="Job configuration file for Galaxy to target.",
        default=None,
        use_global_config=True,
    )


def mulled_containers_option():
    return planemo_option(
        "mulled_containers",
        "--mulled_containers",
        "--biocontainers",
        is_flag=True,
        help="Test tools against mulled containers (forces --docker).",
    )


def install_galaxy_option():
    return planemo_option(
        "--install_galaxy",
        is_flag=True,
        help="Download and configure a disposable copy of Galaxy from github."
    )


def docker_galaxy_image_option():
    return planemo_option(
        "--docker_galaxy_image",
        default="bgruening/galaxy-stable",
        use_global_config=True,
        help=("Docker image identifier for docker-galaxy-flavor used if "
              "engine type is specified as ``docker-galaxy``. Defaults to "
              "to bgruening/galaxy-stable.")
    )


def no_cache_galaxy_option():
    return planemo_option(
        "--no_cache_galaxy",
        is_flag=True,
        help=("Skip caching of Galaxy source and dependencies obtained with "
              "--install_galaxy. Not caching this results in faster "
              "downloads (no git) - so is better on throw away instances such "
              "with TravisCI. ")
    )


def galaxy_branch_option():
    return planemo_option(
        "--galaxy_branch",
        default=None,
        use_global_config=True,
        help=("Branch of Galaxy to target (defaults to master) if a Galaxy "
              "root isn't specified.")
    )


def galaxy_source_option():
    return planemo_option(
        "--galaxy_source",
        default=None,
        use_global_config=True,
        help=("Git source of Galaxy to target (defaults to the official "
              "galaxyproject github source if a Galaxy root isn't "
              "specified.")
    )


def skip_install_option():
    return planemo_option(
        "--skip_install",
        is_flag=True,
        help="Skip installation - only source requirements already available."
    )


def brew_option():
    return planemo_option(
        "--brew",
        use_global_config=True,
        type=click.Path(exists=True, file_okay=True, dir_okay=False),
        help="Homebrew 'brew' executable to use."
    )


def conda_prefix_option():
    return planemo_option(
        "--conda_prefix",
        use_global_config=True,
        use_env_var=True,
        type=click.Path(file_okay=False, dir_okay=True),
        help="Conda prefix to use for conda dependency commands."
    )


def conda_exec_option():
    return planemo_option(
        "--conda_exec",
        use_global_config=True,
        type=click.Path(exists=True, file_okay=True, dir_okay=False),
        help="Location of conda executable."
    )


def conda_debug_option():
    return planemo_option(
        "--conda_debug",
        is_flag=True,
        help="Enable more verbose conda logging."
    )


def conda_use_local_option():
    return planemo_option(
        "--conda_use_local",
        is_flag=True,
        help="Use locally built packages while building Conda environments."
    )


def conda_ensure_channels_option():
    return planemo_option(
        "conda_ensure_channels",
        "--conda_channels",
        "--conda_ensure_channels",
        type=str,
        use_global_config=True,
        use_env_var=True,
        help=("Ensure conda is configured with specified comma separated "
              "list of channels."),
        default="iuc,bioconda,r,defaults,conda-forge",
    )


def conda_copy_dependencies_option():
    return planemo_option(
        "--conda_copy_dependencies",
        is_flag=True,
        help=("Conda dependency resolution for Galaxy will copy dependencies "
              "instead of attempting to link them.")
    )


def conda_auto_install_option():
    return planemo_option(
        "--conda_auto_install/--no_conda_auto_install",
        is_flag=True,
        default=True,
        help=("Conda dependency resolution for Galaxy will auto install "
              "will attempt to install requested but missing packages.")
    )


def conda_auto_init_option():
    return planemo_option(
        "--conda_auto_init/--no_conda_auto_init",
        is_flag=True,
        default=True,
        help=("Conda dependency resolution for Galaxy will auto install "
              "conda itself using miniconda if not availabe on conda_prefix.")
    )


def conda_global_option():
    return planemo_option(
        "--global",
        is_flag=True,
        default=False,
        help=("Install Conda dependencies globally instead of in requirement specific "
              "environments packaged for tools. If the Conda bin directory is on your "
              "PATH, tools may still use binaries but this is more designed for "
              "interactive testing and debugging.")
    )


def required_tool_arg(allow_uris=False):
    """ Decorate click method as requiring the path to a single tool.
    """
    arg_type_class = click.Path if not allow_uris else UriLike
    arg_type = arg_type_class(
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
    )
    if allow_uris:
        name = "uri"
        metavar = "TOOL_URI"
    else:
        name = "path"
        metavar = "TOOL_PATH"
    return click.argument(name, metavar=metavar, type=arg_type)


def required_job_arg():
    """ Decorate click method as requiring the path to a single tool.
    """
    arg_type = click.Path(
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=False,
    )
    return click.argument("job_path", metavar="JOB_PATH", type=arg_type)


def _optional_tools_default(ctx, param, value):
    if param.name in ["paths", "uris"] and len(value) == 0:
        return [os.path.abspath(os.getcwd())]
    else:
        return value


def optional_tools_or_packages_arg(multiple=False):
    """ Decorate click method as optionally taking in the path to a tool
    or directory of tools or a Conda package. If no such argument is given
    the current working directory will be treated as a directory of tools.
    """
    name = "paths" if multiple else "path"
    nargs = -1 if multiple else 1
    return click.argument(
        name,
        metavar="TARGET",
        nargs=nargs,
    )


class UriLike(click.Path):

    def convert(self, value, param, ctx):
        if "://" in value:
            return value
        else:
            return super(UriLike, self).convert(value, param, ctx)


def optional_tools_arg(multiple=False, allow_uris=False):
    """ Decorate click method as optionally taking in the path to a tool
    or directory of tools. If no such argument is given the current working
    directory will be treated as a directory of tools.
    """
    arg_type_class = click.Path if not allow_uris else UriLike
    arg_type = arg_type_class(
        exists=True,
        file_okay=True,
        dir_okay=True,
        readable=True,
        resolve_path=True,
    )
    if allow_uris:
        name = "uris" if multiple else "uri"
    else:
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


def recipe_arg(multiple=True):
    name = "paths" if multiple else "path"
    nargs = -1 if multiple else 1
    return click.argument(
        name,
        metavar="RECIPE_DIR",
        type=click.Path(
            exists=True,
            file_okay=True,
            dir_okay=True,
            resolve_path=True,
        ),
        nargs=nargs,
        callback=_optional_tools_default,
    )


def optional_project_arg(exists=True, default="."):
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
        default=default,
        type=arg_type
    )


def no_cleanup_option():
    return planemo_option(
        "--no_cleanup",
        is_flag=True,
        help=("Do not cleanup temp files created for and by Galaxy.")
    )


def docker_enable_option():
    return planemo_option(
        "--docker/--no_docker",
        default=False,
        help=("Run Galaxy tools in Docker if enabled.")
    )


def docker_cmd_option():
    return planemo_option(
        "--docker_cmd",
        default=docker_util.DEFAULT_DOCKER_COMMAND,
        help="Command used to launch docker (defaults to docker)."
    )


def docker_sudo_option():
    return planemo_option(
        "--docker_sudo/--no_docker_sudo",
        is_flag=True,
        help="Flag to use sudo when running docker."
    )


def docker_sudo_cmd_option():
    return planemo_option(
        "--docker_sudo_cmd",
        help="sudo command to use when --docker_sudo is enabled " +
             "(defaults to sudo).",
        default=docker_util.DEFAULT_SUDO_COMMAND,
        use_global_config=True,
    )


def docker_host_option():
    return planemo_option(
        "--docker_host",
        help="Docker host to target when executing docker commands " +
             "(defaults to localhost).",
        use_global_config=True,
        default=docker_util.DEFAULT_HOST,
    )


def docker_config_options():
    return _compose(
        docker_cmd_option(),
        docker_sudo_option(),
        docker_host_option(),
        docker_sudo_cmd_option(),
    )


def galaxy_docker_options():
    return _compose(
        docker_enable_option(),
        docker_config_options(),
    )


def shed_owner_option():
    return planemo_option(
        "--owner",
        help="Tool Shed repository owner (username)."
    )


def shed_name_option():
    return planemo_option(
        "--name",
        help="Tool Shed repository name (defaults to the inferred "
             "tool directory name)."
    )


def validate_shed_target_callback(ctx, param, value):
    if value is None:
        ctx.fail("default_shed_target set to None, must specify a value for "
                 "--shed_target to run this command.")
    return value


def shed_target_option():
    return planemo_option(
        "-t",
        "--shed_target",
        help="Tool Shed to target (this can be 'toolshed', 'testtoolshed', "
             "'local' (alias for http://localhost:9009/), an arbitrary url "
             "or mappings defined ~/.planemo.yml.",
        default=None,
        use_global_config=True,
        callback=validate_shed_target_callback,
    )


def shed_key_option():
    return planemo_option(
        "--shed_key",
        help=("API key for Tool Shed access. An API key is required unless "
              "e-mail and password is specified. This key can be specified "
              "with either --shed_key or --shed_key_from_env.")
    )


def shed_key_from_env_option():
    return planemo_option(
        "--shed_key_from_env",
        help="Environment variable to read API key for Tool Shed access from."
    )


def shed_email_option():
    return planemo_option(
        "--shed_email",
        help="E-mail for Tool Shed auth (required unless shed_key is "
             "specified)."
    )


def shed_password_option():
    return planemo_option(
        "--shed_password",
        help="Password for Tool Shed auth (required unless shed_key is "
             "specified)."
    )


def shed_skip_upload():
    return planemo_option(
        "--skip_upload",
        is_flag=True,
        help=("Skip upload contents as part of operation, only update "
              "metadata.")
    )


def shed_skip_metadata():
    return planemo_option(
        "--skip_metadata",
        is_flag=True,
        help=("Skip metadata update as part of operation, only upload "
              "new contents.")
    )


def shed_message_option():
    return planemo_option(
        "-m",
        "--message",
        help="Commit message for tool shed upload."
    )


def shed_force_create_option():
    return planemo_option(
        "--force_repository_creation",
        help=("If a repository cannot be found for the specified user/repo "
              "name pair, then automatically create the repository in the "
              "toolshed."),
        is_flag=True,
        default=False
    )


def shed_check_diff_option():
    return planemo_option(
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


def conda_target_options(include_local=True):
    return _compose(
        conda_prefix_option(),
        conda_exec_option(),
        conda_debug_option(),
        conda_ensure_channels_option(),
        conda_use_local_option(),
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
        brew_dependency_resolution(),
        shed_dependency_resolution(),
        no_dependency_resolution(),
        conda_target_options(),
        conda_dependency_resolution(),
        conda_copy_dependencies_option(),
        conda_auto_install_option(),
        conda_auto_init_option(),
        # Profile options...
        profile_option(),
        profile_database_options(),
        file_path_option(),
        database_connection_option(),
        shed_tools_conf_option(),
        shed_tools_directory_option(),
    )


def galaxy_target_options():
    return _compose(
        galaxy_root_option(),
        galaxy_database_seed_option(),
        extra_tools_option(),
        install_galaxy_option(),
        galaxy_branch_option(),
        galaxy_source_option(),
        skip_venv_option(),
        no_cache_galaxy_option(),
        no_cleanup_option(),
        galaxy_email_option(),
        galaxy_docker_options(),
        mulled_containers_option(),
        # Profile options...
        job_config_option(),
        tool_dependency_dir_option(),
    )


def pid_file_option():
    return planemo_option(
        "--pid_file",
        default=None,
        help="Location of pid file is executed with --daemon."
    )


def daemon_option():
    return planemo_option(
        "--daemon",
        is_flag=True,
        help="Serve Galaxy process as a daemon."
    )


def profile_option():
    return planemo_option(
        "--profile",
        type=str,
        default=None,
        help="Location of pid file is executed with --daemon."
    )


def galaxy_serve_options():
    return _compose(
        galaxy_run_options(),
        serve_engine_option(),
        non_strict_cwl_option(),
        docker_galaxy_image_option(),
        galaxy_config_options(),
        daemon_option(),
        pid_file_option(),
    )


def shed_fail_fast_option():
    return planemo_option(
        "--fail_fast",
        is_flag=True,
        default=False,
        help="If multiple repositories are specified and an error occurs "
             "stop immediately instead of processing remaining repositories."
    )


def lint_xsd_option():
    return planemo_option(
        "--xsd/--no_xsd",
        is_flag=True,
        default=True,
        help=("Include tool XSD validation in linting process.")
    )


def report_level_option():
    return planemo_option(
        "--report_level",
        type=click.Choice(["all", "warn", "error"]),
        default="all",
    )


def report_xunit():
    return planemo_option(
        "--report_xunit",
        type=click.Path(file_okay=True, resolve_path=True),
        help="Output an XUnit report, useful for CI testing",
        default=None,
    )


def skip_option():
    return planemo_option(
        "-s",
        "--skip",
        default=None,
        help=("Comma-separated list of lint tests to skip (e.g. passing "
              "--skip 'citations,xml_order' would skip linting of citations "
              "and best-practice XML ordering.")
    )


def fail_level_option():
    return planemo_option(
        "--fail_level",
        type=click.Choice(['warn', 'error']),
        default="warn"
    )


def recursive_shed_option():
    return recursive_option(
        "Recursively perform command for nested repository directories.",
    )


def recursive_option(help="Recursively perform command for subdirectories."):
    return planemo_option(
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


def engine_options():
    return _compose(
        run_engine_option(),
        non_strict_cwl_option(),
        cwltool_no_container_option(),
        docker_galaxy_image_option(),
    )


def test_report_options():
    return _compose(
        planemo_option(
            "--test_output",
            type=click.Path(file_okay=True, resolve_path=True),
            use_global_config=True,
            default="tool_test_output.html",
            help=("Output test report (HTML - for humans) defaults to "
                  "tool_test_output.html."),
        ),
        planemo_option(
            "--test_output_text",
            type=click.Path(file_okay=True, resolve_path=True),
            use_global_config=True,
            help=("Output test report (Basic text - for display in CI)"),
            default=None,
        ),
        planemo_option(
            "--test_output_markdown",
            type=click.Path(file_okay=True, resolve_path=True),
            use_global_config=True,
            help=("Output test report (Markdown style - for humans & "
                  "computers)"),
            default=None,
        ),
    )


def profile_name_argument():
    return click.argument(
        'profile_name',
        metavar="PROFILE_NAME",
        type=str,
    )


def database_identifier_argument():
    return click.argument(
        'identifier',
        metavar="IDENTIFIER",
        type=str,
    )


def postgres_datatype_type_option():
    return planemo_option(
        "--postgres",
        "database_type",
        flag_value="postgres",
        help=("Use postgres database type."),
    )


def database_type_option():
    return planemo_option(
        "--database_type",
        default="postgres",
        type=click.Choice([
            "postgres",
            "sqlite",
        ]),
        use_global_config=True,
        help=("Type of database to use for profile - "
              "currently only 'postgres' is available."),
    )


def database_source_options():
    """Database connection options for commands that utilize a database."""
    return _compose(
        planemo_option(
            "--postgres_psql_path",
            default="psql",
            use_global_config=True,
            help=("Name or or path to postgres client binary (psql)."),
        ),
        planemo_option(
            "--postgres_database_user",
            default="postgres",
            use_global_config=True,
            help=("Postgres username for managed development databases."),
        ),
        planemo_option(
            "--postgres_database_host",
            default=None,
            use_global_config=True,
            help=("Postgres host name for managed development databases."),
        ),
        planemo_option(
            "--postgres_database_port",
            default=None,
            use_global_config=True,
            help=("Postgres port for managed development databases."),
        ),
    )


def profile_database_options():
    return _compose(
        postgres_datatype_type_option(),
        database_type_option(),
        database_source_options(),
    )


def test_options():
    return _compose(
        planemo_option(
            "--update_test_data",
            is_flag=True,
            help="Update test-data directory with job outputs (normally"
                 " written to directory --job_output_files if specified.)"
        ),
        test_report_options(),
        planemo_option(
            "--test_output_xunit",
            type=click.Path(file_okay=True, resolve_path=True),
            use_global_config=True,
            help="Output test report (xUnit style - for computers).",
        ),
        planemo_option(
            "--test_output_json",
            type=click.Path(file_okay=True, resolve_path=True),
            use_global_config=True,
            help=("Output test report (planemo json) defaults to "
                  "tool_test_output.json."),
            default="tool_test_output.json",
        ),
        planemo_option(
            "--job_output_files",
            type=click.Path(file_okay=False, resolve_path=True),
            help="Write job outputs to specified directory.",
            default=None,
        ),
        planemo_option(
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
        planemo_option(
            "--download_cache",
            type=click.Path(file_okay=False, resolve_path=True),
            use_global_config=True,
            help=("Directory to cache downloaded files, default is $DOWNLOAD_CACHE"),
            default=None,
        ),
    )


def filter_exclude_option():
    return planemo_option(
        "--exclude",
        type=click.Path(resolve_path=False),
        multiple=True,
        help="Paths to exclude.",
    )


def filter_exclude_from_option():
    return planemo_option(
        "--exclude_from",
        type=click.Path(exists=True, file_okay=True, dir_okay=False, resolve_path=True),
        multiple=True,
        help="File of paths to exclude.",
    )


def filter_changed_in_commit_option():
    return planemo_option(
        "--changed_in_commit_range",
        help="Exclude paths unchanged in git commit range.",
    )


def ci_chunk_count_option():
    return planemo_option(
        "--chunk_count",
        type=int,
        help="Split output into chunks of this many item and print --chunk such group.",
        default=1,
    )


def ci_chunk_option():
    return planemo_option(
        "--chunk",
        type=int,
        help=("When output is split into --chunk_count groups, output the group 0-indexed"
              "by this option."),
        default=0,
    )


def ci_output_option():
    return planemo_option(
        "--output",
        help="File to output to, or - for standard output.",
        default="-",
    )


def ci_find_options():
    return _compose(
        filter_exclude_option(),
        filter_exclude_from_option(),
        filter_changed_in_commit_option(),
        ci_chunk_count_option(),
        ci_chunk_option(),
        ci_output_option(),
    )


def tool_init_id_option(prompt=True):
    return planemo_option(
        "-i",
        "--id",
        type=click.STRING,
        prompt=prompt,
        help="Short identifier for new tool (no whitespace)",
    )


def tool_init_tool_option():
    return planemo_option(
        "-t",
        "--tool",
        default=None,
        type=click.Path(exists=False,
                        file_okay=True,
                        dir_okay=False,
                        writable=True,
                        resolve_path=True),
        help="Output path for new tool (default is <id>.xml)",
    )


def tool_init_name_option(prompt=True, help="Name for new tool (user facing)"):
    return planemo_option(
        "-n",
        "--name",
        type=click.STRING,
        prompt=prompt,
        help=help,
    )


def tool_init_version_option():
    return planemo_option(
        "--version",
        default="0.1.0",
        type=click.STRING,
        help="Tool XML version.",
    )


def tool_init_description_option():
    return planemo_option(
        "-d",
        "--description",
        type=click.STRING,
        default=None,
        prompt=False,
        help="Short description for new tool (user facing)",
    )


def tool_init_command_option():
    return planemo_option(
        "-c",
        "--command",
        type=click.STRING,
        default=None,
        prompt=False,
        help=("Command potentially including cheetah variables ()"
              "(e.g. 'seqtk seq -a $input > $output')"),
    )


def tool_init_doi_option():
    return planemo_option(
        "--doi",
        type=click.STRING,
        default=None,
        multiple=True,
        prompt=False,
        help=("Supply a DOI (http://www.doi.org/) easing citation of the tool "
              "for Galxy users (e.g. 10.1101/014043).")
    )


def tool_init_test_case_option():
    return planemo_option(
        "--test_case",
        is_flag=True,
        default=None,
        prompt=False,
        help=("For use with --example_commmand, generate a tool test case from "
              "the supplied example."),
    )


def tool_init_macros_option():
    return planemo_option(
        "--macros",
        is_flag=True,
        default=None,
        prompt=False,
        help="Generate a macros.xml for reuse across many tools.",
    )


def tool_init_cite_url_option():
    return planemo_option(
        "--cite_url",
        type=click.STRING,
        default=None,
        multiple=True,
        prompt=False,
        help=("Supply a URL for citation.")
    )


def tool_init_input_option():
    return planemo_option(
        "--input",
        type=click.STRING,
        default=None,
        prompt=False,
        multiple=True,
        help="An input description (e.g. input.fasta)",
    )


def tool_init_output_option():
    return planemo_option(
        "--output",
        type=click.STRING,
        multiple=True,
        default=None,
        prompt=False,
        help=("An output location (e.g. output.bam), the Galaxy datatype is "
              "inferred from the extension."),
    )


def tool_init_help_text_option():
    return planemo_option(
        "--help_text",
        type=click.STRING,
        default=None,
        prompt=False,
        help="Help text (reStructuredText)",
    )


def tool_init_help_from_command_option():
    return planemo_option(
        "--help_from_command",
        type=click.STRING,
        default=None,
        prompt=False,
        help="Auto populate help from supplied command.",
    )


def tool_init_example_input_option():
    return planemo_option(
        "--example_input",
        type=click.STRING,
        default=None,
        prompt=False,
        multiple=True,
        help=("For use with --example_command, replace input file (e.g. 2.fastq "
              "with a data input parameter)."),
    )


def tool_init_example_output_option():
    return planemo_option(
        "--example_output",
        type=click.STRING,
        default=None,
        prompt=False,
        multiple=True,
        help=("For use with --example_command, replace input file (e.g. 2.fastq "
              "with a tool output)."),
    )


def tool_init_named_output_option():
    return planemo_option(
        "--named_output",
        type=click.STRING,
        multiple=True,
        default=None,
        prompt=False,
        help=("Create a named output for use with command block for example "
              "specify --named_output=output1.bam and then use '-o $output1' "
              "in your command block."),
    )


def tool_init_version_command_option():
    return planemo_option(
        "--version_command",
        type=click.STRING,
        default=None,
        prompt=False,
        help="Command to print version (e.g. 'seqtk --version')",
    )


REQUIREMENT_HELP = "Add a tool requirement package (e.g. 'seqtk' or 'seqtk@1.68')."


def tool_init_requirement_option(help=REQUIREMENT_HELP):
    return planemo_option(
        "--requirement",
        type=click.STRING,
        default=None,
        multiple=True,
        prompt=False,
        help=help,
    )


def tool_init_container_option():
    return planemo_option(
        "--container",
        type=click.STRING,
        default=None,
        multiple=True,
        prompt=False,
        help="Add a Docker image identifier for this tool."
    )


EXAMPLE_COMMAND_HELP = (
    "Example to command with paths to build Cheetah template from "
    "(e.g. 'seqtk seq -a 2.fastq > 2.fasta'). Option cannot be used "
    "with --command, should be used --example_input and "
    "--example_output."
)


def tool_init_example_command_option(help=EXAMPLE_COMMAND_HELP):
    return planemo_option(
        "--example_command",
        type=click.STRING,
        default=None,
        prompt=False,
        help=help,
    )


def mulled_conda_option():
    return planemo_option(
        "--mulled_conda_version",
        type=click.STRING,
        default=None,
        help=("Install a specific version of Conda before running the command, by "
              "default the version that comes with the continuumio miniconda3 image "
              "will be used under Linux and under Mac OS X Conda will be upgraded to "
              "to work around a bug in 4.2.")
    )


def mulled_namespace_option():
    return planemo_option(
        "--mulled_namespace",
        type=click.STRING,
        default="biocontainers",
        help=("Build a mulled image with the specified namespace - defaults to "
              "biocontainers. Galaxy currently only recognizes images with the "
              "namespace biocontainers.")
    )


def mulled_action_option():
    return planemo_option(
        "--mulled_command",
        type=click.STRING,
        default="build",
        help=("Mulled action to perform for targets - this defaults to 'build'. "
              "Set this to build-and-test to also test the resulting container.")
    )


def mulled_options():
    return _compose(
        mulled_conda_option(),
        mulled_namespace_option(),
        mulled_action_option(),
    )
