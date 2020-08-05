"""Abstractions for setting up a Galaxy instance."""
from __future__ import absolute_import
from __future__ import print_function

import abc
import contextlib
import os
import random
import shutil
from string import Template
from tempfile import mkdtemp

from galaxy.containers.docker_model import DockerVolume
from galaxy.tool_util.deps import docker_util
from galaxy.util.commands import argv_to_str
from pkg_resources import parse_version
from six import (
    add_metaclass,
    iteritems
)
from six.moves import shlex_quote

from planemo import git
from planemo.config import OptionSource
from planemo.deps import ensure_dependency_resolvers_conf_configured
from planemo.docker import docker_host_args
from planemo.io import (
    communicate,
    kill_pid_file,
    shell,
    shell_join,
    untar_to,
    wait_on,
    warn,
    write_file,
)
from planemo.mulled import build_involucro_context
from planemo.shed import tool_shed_url
from planemo.virtualenv import DEFAULT_PYTHON_VERSION
from .api import (
    DEFAULT_MASTER_API_KEY,
    gi,
    user_api_key,
)
from .distro_tools import (
    DISTRO_TOOLS_ID_TO_PATH
)
from .run import (
    setup_common_startup_args,
    setup_venv,
)
from .workflows import (
    find_tool_ids,
    import_workflow,
    install_shed_repos,
)


NO_TEST_DATA_MESSAGE = (
    "planemo couldn't find a target test-data directory, you should likely "
    "create a test-data directory or pass an explicit path using --test_data."
)

WEB_SERVER_CONFIG_TEMPLATE = """
[server:${server_name}]
use = egg:Paste#http
port = ${port}
host = ${host}
use_threadpool = True
threadpool_kill_thread_limit = 10800
[app:main]
paste.app_factory = galaxy.web.buildapp:app_factory
static_dir = static/
"""

TOOL_CONF_TEMPLATE = """<toolbox>
  <tool file="data_source/upload.xml" />
  ${tool_definition}
</toolbox>
"""

SHED_TOOL_CONF_TEMPLATE = """<?xml version="1.0"?>
<toolbox tool_path="${shed_tool_path}">
</toolbox>
"""

SHED_DATA_MANAGER_CONF_TEMPLATE = """<?xml version="1.0"?>
<data_managers>
</data_managers>
"""

EMPTY_JOB_METRICS_TEMPLATE = """<?xml version="1.0"?>
<job_metrics>
</job_metrics>
"""

TOOL_SHEDS_CONF = """<tool_sheds>
  <tool_shed name="Target Shed" url="${shed_target_url}" />
</tool_sheds>
"""

JOB_CONFIG_LOCAL = """<job_conf>
    <plugins>
        <plugin id="planemo_runner" type="runner" load="galaxy.jobs.runners.local:LocalJobRunner" workers="4"/>
    </plugins>
    <handlers>
        <handler id="main"/>
    </handlers>
    <destinations default="planemo_dest">
        <destination id="planemo_dest" runner="planemo_runner">
            <param id="require_container">${require_container}</param>
            <param id="docker_enabled">${docker_enable}</param>
            <param id="docker_sudo">${docker_sudo}</param>
            <param id="docker_sudo_cmd">${docker_sudo_cmd}</param>
            <param id="docker_cmd">${docker_cmd}</param>
            ${docker_host_param}
        </destination>
        <destination id="upload_dest" runner="planemo_runner">
            <param id="docker_enable">false</param>
        </destination>
    </destinations>
    <tools>
        <tool id="upload1" destination="upload_dest" />
    </tools>
</job_conf>
"""

LOGGING_TEMPLATE = """
## Configure Python loggers.
[loggers]
keys = root,paste,displayapperrors,galaxydeps,galaxymasterapikey,galaxy

[handlers]
keys = console

[formatters]
keys = generic

[logger_root]
level = WARN
handlers = console

[logger_paste]
level = WARN
handlers = console
qualname = paste
propagate = 0

[logger_galaxydeps]
level = DEBUG
handlers = console
qualname = galaxy.tools.deps
propagate = 0

[logger_galaxymasterapikey]
level = WARN
handlers = console
qualname = galaxy.web.framework.webapp
propagate = 0

[logger_displayapperrors]
level = ERROR
handlers =
qualname = galaxy.datatypes.display_applications.application
propagate = 0

[logger_galaxy]
level = ${log_level}
handlers = console
qualname = galaxy
propagate = 0

[handler_console]
class = StreamHandler
args = (sys.stderr,)
level = DEBUG
formatter = generic

[formatter_generic]
format = %(asctime)s %(levelname)-5.5s [%(name)s] %(message)s
"""

REFGENIE_CONFIG_TEMPLATE = """
config_version: 0.3
genome_folder: '%s'
genome_servers: ['http://refgenomes.databio.org']
genomes: null
"""

EMPTY_TOOL_CONF_TEMPLATE = """<toolbox></toolbox>"""

DEFAULT_GALAXY_BRANCH = "master"
DEFAULT_GALAXY_SOURCE = "https://github.com/galaxyproject/galaxy"
CWL_GALAXY_SOURCE = "https://github.com/common-workflow-language/galaxy"

DATABASE_LOCATION_TEMPLATE = "sqlite:///%s?isolation_level=IMMEDIATE"

COMMAND_STARTUP_COMMAND = './scripts/common_startup.sh ${COMMON_STARTUP_ARGS}'

CLEANUP_IGNORE_ERRORS = True
DEFAULT_GALAXY_BRAND = 'Configured by Planemo'
DEFAULT_TOOL_INSTALL_TIMEOUT = 60 * 60 * 1
UNINITIALIZED = object()


@contextlib.contextmanager
def galaxy_config(ctx, runnables, **kwds):
    """Set up a ``GalaxyConfig`` in an auto-cleaned context."""
    c = local_galaxy_config
    if kwds.get("dockerize", False):
        c = docker_galaxy_config
    elif kwds.get("external", False):
        c = external_galaxy_config

    with c(ctx, runnables, **kwds) as config:
        yield config


def simple_docker_volume(path):
    path = os.path.abspath(path)
    return DockerVolume("%s:%s:rw" % (path, path))


@contextlib.contextmanager
def docker_galaxy_config(ctx, runnables, for_tests=False, **kwds):
    """Set up a ``GalaxyConfig`` for Docker container."""
    test_data_dir = _find_test_data(runnables, **kwds)

    with _config_directory(ctx, **kwds) as config_directory:
        def config_join(*args):
            return os.path.join(config_directory, *args)

        ensure_dependency_resolvers_conf_configured(ctx, kwds, os.path.join(config_directory, "resolvers_conf.xml"))
        _handle_job_metrics(config_directory, kwds)
        _handle_refgenie_config(config_directory, kwds)

        shed_tool_conf = "config/shed_tool_conf.xml"
        all_tool_paths = _all_tool_paths(runnables, **kwds)

        tool_directories = set([])  # Things to mount...
        for tool_path in all_tool_paths:
            directory = os.path.dirname(os.path.normpath(tool_path))
            if os.path.exists(directory):
                tool_directories.add(directory)

        # TODO: remap these.
        tool_volumes = []
        for tool_directory in tool_directories:
            volume = simple_docker_volume(tool_directory)
            tool_volumes.append(volume)

        empty_tool_conf = config_join("empty_tool_conf.xml")

        tool_conf = config_join("tool_conf.xml")

        shed_tool_path = kwds.get("shed_tool_path") or config_join("shed_tools")
        _ensure_directory(shed_tool_path)

        sheds_config_path = _configure_sheds_config_file(
            ctx, config_directory, **kwds
        )
        port = _get_port(kwds)
        properties = _shared_galaxy_properties(config_directory, kwds, for_tests=for_tests)
        _handle_container_resolution(ctx, kwds, properties)
        master_api_key = _get_master_api_key(kwds)

        template_args = dict(
            shed_tool_path=shed_tool_path,
            tool_conf=tool_conf,
        )
        tool_config_file = "%s,%s" % (tool_conf, shed_tool_conf)

        _write_tool_conf(ctx, all_tool_paths, tool_conf)
        write_file(empty_tool_conf, EMPTY_TOOL_CONF_TEMPLATE)

        properties.update(dict(
            tool_config_file=tool_config_file,
            tool_sheds_config_file=sheds_config_path,
            migrated_tools_config=empty_tool_conf,
        ))

        server_name = "planemo%d" % random.randint(0, 100000)

        # Value substitutions in Galaxy properties - for consistency with
        # non-Dockerized version.
        template_args = dict(
        )
        env = _build_env_for_galaxy(properties, template_args)
        env["NONUSE"] = "nodejs,proftp,reports"
        if ctx.verbose:
            env["GALAXY_LOGGING"] = "full"

        # TODO: setup FTP upload dir and disable FTP server in container.
        _build_test_env(properties, env)

        docker_target_kwds = docker_host_args(**kwds)
        volumes = tool_volumes + [simple_docker_volume(config_directory)]
        export_directory = kwds.get("export_directory", None)
        if export_directory is not None:
            volumes.append(DockerVolume("%s:/export:rw" % export_directory))

        # TODO: Allow this to real Docker volumes and allow multiple.
        extra_volume = kwds.get("docker_extra_volume")
        if extra_volume:
            volumes.append(simple_docker_volume(extra_volume))
        yield DockerGalaxyConfig(
            ctx,
            config_directory,
            env,
            test_data_dir,
            port,
            server_name,
            master_api_key,
            runnables,
            docker_target_kwds=docker_target_kwds,
            volumes=volumes,
            export_directory=export_directory,
            kwds=kwds,
        )


@contextlib.contextmanager
def local_galaxy_config(ctx, runnables, for_tests=False, **kwds):
    """Set up a ``GalaxyConfig`` in an auto-cleaned context."""
    test_data_dir = _find_test_data(runnables, **kwds)
    tool_data_table = _find_tool_data_table(
        runnables,
        test_data_dir=test_data_dir,
        **kwds
    )
    data_manager_config_paths = [r.data_manager_conf_path for r in runnables if r.data_manager_conf_path]
    galaxy_root = _find_galaxy_root(ctx, **kwds)
    install_galaxy = kwds.get("install_galaxy", False)
    if galaxy_root is not None:
        if os.path.isdir(galaxy_root) and not os.listdir(galaxy_root):
            os.rmdir(galaxy_root)
        if os.path.isdir(galaxy_root) and install_galaxy:
            raise Exception("%s is an existing non-empty directory, cannot install Galaxy again" % galaxy_root)

    # Duplicate block in docker variant above.
    if kwds.get("mulled_containers", False) and not kwds.get("docker", False):
        if ctx.get_option_source("docker") != OptionSource.cli:
            kwds["docker"] = True
        else:
            raise Exception("Specified no docker and mulled containers together.")

    with _config_directory(ctx, **kwds) as config_directory:
        def config_join(*args):
            return os.path.join(config_directory, *args)

        install_env = {}
        if kwds.get('galaxy_skip_client_build', True):
            install_env['GALAXY_SKIP_CLIENT_BUILD'] = '1'
        if galaxy_root is None:
            galaxy_root = config_join("galaxy-dev")
        if not os.path.isdir(galaxy_root):
            _build_eggs_cache(ctx, install_env, kwds)
            _install_galaxy(ctx, galaxy_root, install_env, kwds)

        if parse_version(kwds.get('galaxy_python_version') or DEFAULT_PYTHON_VERSION) >= parse_version('3'):
            # on python 3 we use gunicorn,
            # which requires 'main' as server name
            server_name = 'main'
        else:
            server_name = "planemo%d" % random.randint(0, 100000)
        # Once we don't have to support earlier than 18.01 - try putting these files
        # somewhere better than with Galaxy.
        log_file = "%s.log" % server_name
        pid_file = "%s.pid" % server_name
        ensure_dependency_resolvers_conf_configured(ctx, kwds, os.path.join(config_directory, "resolvers_conf.xml"))
        _handle_job_config_file(config_directory, server_name, kwds)
        _handle_job_metrics(config_directory, kwds)
        _handle_refgenie_config(config_directory, kwds)
        file_path = kwds.get("file_path") or config_join("files")
        _ensure_directory(file_path)

        tool_dependency_dir = kwds.get("tool_dependency_dir") or config_join("deps")
        _ensure_directory(tool_dependency_dir)

        shed_tool_conf = kwds.get("shed_tool_conf") or config_join("shed_tools_conf.xml")
        all_tool_paths = _all_tool_paths(runnables, **kwds)
        empty_tool_conf = config_join("empty_tool_conf.xml")

        tool_conf = config_join("tool_conf.xml")

        shed_data_manager_config_file = config_join("shed_data_manager_conf.xml")

        shed_tool_path = kwds.get("shed_tool_path") or config_join("shed_tools")
        _ensure_directory(shed_tool_path)

        sheds_config_path = _configure_sheds_config_file(
            ctx, config_directory, **kwds
        )

        database_location = config_join("galaxy.sqlite")
        master_api_key = _get_master_api_key(kwds)
        dependency_dir = os.path.join(config_directory, "deps")
        _ensure_directory(shed_tool_path)
        port = _get_port(kwds)
        template_args = dict(
            port=port,
            host=kwds.get("host", "127.0.0.1"),
            server_name=server_name,
            temp_directory=config_directory,
            shed_tool_path=shed_tool_path,
            database_location=database_location,
            tool_conf=tool_conf,
            debug=kwds.get("debug", "true"),
            id_secret=kwds.get("id_secret", "test_secret"),
            log_level="DEBUG" if ctx.verbose else "INFO",
        )
        tool_config_file = "%s,%s" % (tool_conf, shed_tool_conf)
        # Setup both galaxy_email and older test user test@bx.psu.edu
        # as admins for command_line, etc...
        properties = _shared_galaxy_properties(config_directory, kwds, for_tests=for_tests)
        properties.update(dict(
            server_name="main",
            ftp_upload_dir_template="${ftp_upload_dir}",
            ftp_upload_purge="False",
            ftp_upload_dir=test_data_dir or os.path.abspath('.'),
            ftp_upload_site="Test Data",
            check_upload_content="False",
            tool_dependency_dir=dependency_dir,
            file_path=file_path,
            new_file_path="${temp_directory}/tmp",
            tool_config_file=tool_config_file,
            tool_sheds_config_file=sheds_config_path,
            manage_dependency_relationships="False",
            job_working_directory="${temp_directory}/job_working_directory",
            template_cache_path="${temp_directory}/compiled_templates",
            citation_cache_type="file",
            citation_cache_data_dir="${temp_directory}/citations/data",
            citation_cache_lock_dir="${temp_directory}/citations/lock",
            database_auto_migrate="True",
            enable_beta_tool_formats="True",
            id_secret="${id_secret}",
            log_level="${log_level}",
            debug="${debug}",
            watch_tools="auto",
            default_job_shell="/bin/bash",  # For conda dependency resolution
            tool_data_table_config_path=tool_data_table,
            data_manager_config_file=",".join(data_manager_config_paths) or None,  # without 'or None' may raise IOError in galaxy (see #946)
            integrated_tool_panel_config=("${temp_directory}/"
                                          "integrated_tool_panel_conf.xml"),
            migrated_tools_config=empty_tool_conf,
            test_data_dir=test_data_dir,  # TODO: make gx respect this
            shed_data_manager_config_file=shed_data_manager_config_file,
        ))
        _handle_container_resolution(ctx, kwds, properties)
        write_file(config_join("logging.ini"), _sub(LOGGING_TEMPLATE, template_args))
        properties["database_connection"] = _database_connection(database_location, **kwds)

        _handle_kwd_overrides(properties, kwds)

        # TODO: consider following property
        # watch_tool = False
        # datatypes_config_file = config/datatypes_conf.xml
        # welcome_url = /static/welcome.html
        # logo_url = /
        # sanitize_all_html = True
        # serve_xss_vulnerable_mimetypes = False
        # track_jobs_in_database = None
        # outputs_to_working_directory = False
        # retry_job_output_collection = 0

        env = _build_env_for_galaxy(properties, template_args)
        env.update(install_env)
        _build_test_env(properties, env)
        env['GALAXY_TEST_SHED_TOOL_CONF'] = shed_tool_conf
        env['GALAXY_TEST_DBURI'] = properties["database_connection"]

        env["GALAXY_TEST_UPLOAD_ASYNC"] = "false"
        env["GALAXY_TEST_LOGGING_CONFIG"] = config_join("logging.ini")
        env["GALAXY_DEVELOPMENT_ENVIRONMENT"] = "1"
        # Following are needed in 18.01 to prevent Galaxy from changing log and pid.
        # https://github.com/galaxyproject/planemo/issues/788
        env["GALAXY_LOG"] = log_file
        env["GALAXY_PID"] = pid_file
        web_config = _sub(WEB_SERVER_CONFIG_TEMPLATE, template_args)
        write_file(config_join("galaxy.ini"), web_config)
        _write_tool_conf(ctx, all_tool_paths, tool_conf)
        write_file(empty_tool_conf, EMPTY_TOOL_CONF_TEMPLATE)

        shed_tool_conf_contents = _sub(SHED_TOOL_CONF_TEMPLATE, template_args)
        # Write a new shed_tool_conf.xml if needed.
        write_file(shed_tool_conf, shed_tool_conf_contents, force=False)

        write_file(shed_data_manager_config_file, SHED_DATA_MANAGER_CONF_TEMPLATE)

        yield LocalGalaxyConfig(
            ctx,
            config_directory,
            env,
            test_data_dir,
            port,
            server_name,
            master_api_key,
            runnables,
            galaxy_root,
            kwds,
        )


def _all_tool_paths(runnables, **kwds):
    tool_paths = [r.path for r in runnables if r.has_tools and not r.data_manager_conf_path]
    all_tool_paths = list(tool_paths) + list(kwds.get("extra_tools", []))
    for runnable in runnables:
        if runnable.type.name == "galaxy_workflow":
            tool_ids = find_tool_ids(runnable.path)
            for tool_id in tool_ids:
                if tool_id in DISTRO_TOOLS_ID_TO_PATH:
                    all_tool_paths.append(DISTRO_TOOLS_ID_TO_PATH[tool_id])

    return all_tool_paths


def _shared_galaxy_properties(config_directory, kwds, for_tests):
    """Setup properties useful for local and Docker Galaxy instances.

    Most things related to paths, etc... are very different between Galaxy
    modalities and many taken care of internally to the container in that mode.
    But this method sets up API stuff, tool, and job stuff that can be shared.
    """
    master_api_key = _get_master_api_key(kwds)
    user_email = _user_email(kwds)
    properties = {
        'master_api_key': master_api_key,
        'admin_users': "%s,test@bx.psu.edu" % user_email,
        'expose_dataset_path': "True",
        'cleanup_job': 'never',
        'collect_outputs_from': "job_working_directory",
        'allow_path_paste': "True",
        'check_migrate_tools': "False",
        'use_cached_dependency_manager': str(kwds.get("conda_auto_install", False)),
        'brand': kwds.get("galaxy_brand", DEFAULT_GALAXY_BRAND),
        'strict_cwl_validation': str(not kwds.get("non_strict_cwl", False)),
    }
    if kwds.get("galaxy_single_user", True):
        properties['single_user'] = user_email

    if for_tests:
        empty_dir = os.path.join(config_directory, "empty")
        _ensure_directory(empty_dir)
        properties["tour_config_dir"] = empty_dir
        properties["interactive_environment_plugins_directory"] = empty_dir
        properties["visualization_plugins_directory"] = empty_dir
        properties["refgenie_config_file"] = kwds.get('refgenie_config_file', '')
    return properties


@contextlib.contextmanager
def external_galaxy_config(ctx, runnables, for_tests=False, **kwds):
    yield BaseGalaxyConfig(
        ctx=ctx,
        galaxy_url=kwds.get("galaxy_url", None),
        master_api_key=_get_master_api_key(kwds),
        user_api_key=kwds.get("galaxy_user_key", None),
        runnables=runnables,
        kwds=kwds
    )


def _get_master_api_key(kwds):
    master_api_key = kwds.get("galaxy_admin_key") or DEFAULT_MASTER_API_KEY
    return master_api_key


def _get_port(kwds):
    port = int(kwds.get("port", 9090))
    return port


def _user_email(kwds):
    user_email = kwds.get("galaxy_email")
    return user_email


@contextlib.contextmanager
def _config_directory(ctx, **kwds):
    config_directory = kwds.get("config_directory", None)
    created_config_directory = False
    if not config_directory:
        created_config_directory = True
        config_directory = os.path.realpath(mkdtemp())
        ctx.vlog("Created directory for Galaxy configuration [%s]" % config_directory)
    try:
        yield config_directory
    finally:
        cleanup = not kwds.get("no_cleanup", False)
        if created_config_directory and cleanup:
            shutil.rmtree(config_directory)


@add_metaclass(abc.ABCMeta)
class GalaxyInterface(object):
    """Abstraction around a Galaxy instance.

    Description of a Galaxy instance and how to interact with it - this could
    potentially be a remote, already running instance or an instance Planemo manages
    to execute some task(s).
    """

    @abc.abstractproperty
    def gi(self):
        """Return an admin bioblend Galaxy instance for API interactions."""

    @abc.abstractproperty
    def user_gi(self):
        """Return a user-backed bioblend Galaxy instance for API interactions."""

    @abc.abstractmethod
    def install_repo(self, *args, **kwds):
        """Install specified tool shed repository."""

    @abc.abstractproperty
    def tool_shed_client(self):
        """Return a admin bioblend tool shed client."""

    @abc.abstractmethod
    def wait_for_all_installed(self):
        """Wait for all queued up repositories installs to complete."""

    @abc.abstractmethod
    def install_workflows(self):
        """Install all workflows configured with these planemo arguments."""

    @abc.abstractmethod
    def workflow_id(self, path):
        """Get installed workflow API ID for input path."""

    @abc.abstractproperty
    def version_major(self):
        """Return target Galaxy version."""

    @abc.abstractproperty
    def user_api_config(self):
        """Return the API indicated configuration for user session.

        Calling .config.get_config() with admin GI session would yield
        a different object (admins have different view of Galaxy's
        configuration).
        """

    @property
    def user_is_admin(self):
        return self.user_api_config["is_admin_user"]


@add_metaclass(abc.ABCMeta)
class GalaxyConfig(GalaxyInterface):
    """Specialization of GalaxyInterface for Galaxy instances Planemo manages itself.

    This assumes more than an API connection is available - Planemo needs to be able to
    start and stop the Galaxy instance, recover logs, etc... There are currently two
    implementations - a locally executed Galaxy and one running inside a Docker containe
    """

    @abc.abstractproperty
    def kill(self):
        """Stop the running instance."""

    @abc.abstractmethod
    def startup_command(self, ctx, **kwds):
        """Return a shell command used to startup this instance.

        Among other common planmo kwds, this should respect the
        ``daemon`` keyword.
        """

    @abc.abstractproperty
    def log_contents(self):
        """Retrieve text of log for running Galaxy instance."""

    @abc.abstractmethod
    def cleanup(self):
        """Cleanup allocated resources to run this instance."""

    @abc.abstractproperty
    def use_path_paste(self):
        """Use path paste to upload data.

        This will only be an option if the target user key is an
        admin user key.
        """


class BaseGalaxyConfig(GalaxyInterface):

    def __init__(
        self,
        ctx,
        galaxy_url,
        master_api_key,
        user_api_key,
        runnables,
        kwds,
    ):
        self._ctx = ctx
        self.galaxy_url = galaxy_url
        self.master_api_key = master_api_key
        self._user_api_key = user_api_key
        self.runnables = runnables
        self._kwds = kwds
        self._workflow_ids = {}

        self._target_version = UNINITIALIZED
        self._target_user_config = UNINITIALIZED

    @property
    def gi(self):
        assert self.galaxy_url
        return gi(url=self.galaxy_url, key=self.master_api_key)

    @property
    def user_gi(self):
        user_api_key = self.user_api_key
        assert user_api_key
        return self._gi_for_key(user_api_key)

    @property
    def user_api_key(self):
        # TODO: thread-safe
        if self._user_api_key is None:
            # TODO: respect --galaxy_email - seems like a real bug
            self._user_api_key = user_api_key(self.gi)

        return self._user_api_key

    def _gi_for_key(self, key):
        assert self.galaxy_url
        return gi(url=self.galaxy_url, key=key)

    def install_repo(self, *args, **kwds):
        self.tool_shed_client.install_repository_revision(
            *args, **kwds
        )

    @property
    def tool_shed_client(self):
        return self.gi.toolShed

    def wait_for_all_installed(self):
        def status_ready(repo):
            status = repo["status"]
            if status in ["Installing", "New"]:
                return None
            if status == "Installed":
                return True
            raise Exception("Error installing repo status is %s" % status)

        def ready():
            repos = self.tool_shed_client.get_repositories()
            ready = all(map(status_ready, repos))
            return ready or None

        wait_on(ready, "galaxy tool installation", timeout=DEFAULT_TOOL_INSTALL_TIMEOUT)

    def install_workflows(self):
        for runnable in self.runnables:
            if runnable.type.name in ["galaxy_workflow", "cwl_workflow"]:
                self._install_workflow(runnable)

    def _install_workflow(self, runnable):
        if self._kwds.get("shed_install"):
            install_shed_repos(runnable,
                               self.gi,
                               self._kwds.get("ignore_dependency_problems", False),
                               self._kwds.get("install_tool_dependencies", False),
                               self._kwds.get("install_resolver_dependencies", True),
                               self._kwds.get("install_repository_dependencies", True))

        default_from_path = self._kwds.get("workflows_from_path", False)
        # TODO: Allow serialization so this doesn't need to assume a
        # shared filesystem with Galaxy server.
        from_path = default_from_path or (runnable.type.name == "cwl_workflow")
        workflow = import_workflow(
            runnable.path, admin_gi=self.gi, user_gi=self.user_gi, from_path=from_path
        )
        self._workflow_ids[runnable.path] = workflow["id"]

    def workflow_id(self, path):
        return self._workflow_ids[path]

    @property
    def use_path_paste(self):
        option = self._kwds.get("paste_test_data_paths")
        if option is None:
            return self.default_use_path_paste
        else:
            return option

    @property
    def default_use_path_paste(self):
        return False

    @property
    def version_major(self):
        """Return target Galaxy version."""
        if self._target_version is UNINITIALIZED:
            self._target_version = self.user_gi.config.get_version()["version_major"]
        return self._target_version

    @property
    def user_api_config(self):
        """Return the API indicated configuration for user session."""
        if self._target_user_config is UNINITIALIZED:
            self._target_user_config = self.user_gi.config.get_config()
        return self._target_user_config


class BaseManagedGalaxyConfig(BaseGalaxyConfig):

    def __init__(
        self,
        ctx,
        config_directory,
        env,
        test_data_dir,
        port,
        server_name,
        master_api_key,
        runnables,
        kwds,
    ):
        galaxy_url = "http://localhost:%d" % port
        super(BaseManagedGalaxyConfig, self).__init__(
            ctx=ctx,
            galaxy_url=galaxy_url,
            master_api_key=master_api_key,
            user_api_key=None,
            runnables=runnables,
            kwds=kwds
        )
        self.config_directory = config_directory
        self.env = env
        self.test_data_dir = test_data_dir
        self.port = port
        self.server_name = server_name


class DockerGalaxyConfig(BaseManagedGalaxyConfig):
    """A :class:`GalaxyConfig` description of a Dockerized Galaxy instance."""

    def __init__(
        self,
        ctx,
        config_directory,
        env,
        test_data_dir,
        port,
        server_name,
        master_api_key,
        runnables,
        docker_target_kwds,
        volumes,
        export_directory,
        kwds,
    ):
        super(DockerGalaxyConfig, self).__init__(
            ctx,
            config_directory,
            env,
            test_data_dir,
            port,
            server_name,
            master_api_key,
            runnables,
            kwds,
        )
        self.docker_target_kwds = docker_target_kwds
        self.volumes = volumes
        self.export_directory = export_directory

    def kill(self):
        """Kill planemo container..."""
        kill_command = docker_util.kill_command(
            self.server_name,
            **self.docker_target_kwds
        )
        return shell(kill_command)

    def startup_command(self, ctx, **kwds):
        """Return a shell command used to startup this instance.

        Among other common planmo kwds, this should respect the
        ``daemon`` keyword.
        """
        daemon = kwds.get("daemon", False)
        daemon_str = "" if not daemon else " -d"
        docker_run_extras = "-p %s:80%s" % (self.port, daemon_str)
        env_directives = ["%s='%s'" % item for item in self.env.items()]
        image = kwds.get("docker_galaxy_image", "bgruening/galaxy-stable")
        run_command = docker_util.build_docker_run_command(
            "", image,
            interactive=False,
            env_directives=env_directives,
            working_directory=None,
            name=self.server_name,
            run_extra_arguments=docker_run_extras,
            set_user=False,
            volumes=self.volumes,
            **self.docker_target_kwds
        )
        chmod_command = [
            "chmod",
            "-R",
            "o+rwx",
            self.config_directory,
        ]
        if self.export_directory:
            chmod_command.append(self.export_directory)

        return shell_join(
            argv_to_str(chmod_command),
            run_command,
        )

    @property
    def log_contents(self):
        logs_command = docker_util.logs_command(
            self.server_name,
            **self.docker_target_kwds
        )
        output, _ = communicate(
            logs_command
        )
        return output

    def cleanup(self):
        shutil.rmtree(self.config_directory, CLEANUP_IGNORE_ERRORS)


class LocalGalaxyConfig(BaseManagedGalaxyConfig):
    """A local, non-containerized implementation of :class:`GalaxyConfig`."""

    def __init__(
        self,
        ctx,
        config_directory,
        env,
        test_data_dir,
        port,
        server_name,
        master_api_key,
        runnables,
        galaxy_root,
        kwds,
    ):
        super(LocalGalaxyConfig, self).__init__(
            ctx,
            config_directory,
            env,
            test_data_dir,
            port,
            server_name,
            master_api_key,
            runnables,
            kwds,
        )
        self.galaxy_root = galaxy_root

    def kill(self):
        if self._ctx.verbose:
            shell(["ps", "ax"])
            exists = os.path.exists(self.pid_file)
            print("Killing pid file [%s]" % self.pid_file)
            print("pid_file exists? [%s]" % exists)
            if exists:
                print("pid_file contents are [%s]" % open(self.pid_file, "r").read())
        kill_pid_file(self.pid_file)

    def startup_command(self, ctx, **kwds):
        """Return a shell command used to startup this instance.

        Among other common planemo kwds, this should respect the
        ``daemon`` keyword.
        """
        daemon = kwds.get("daemon", False)
        # TODO: Allow running dockerized Galaxy here instead.
        setup_venv_command = setup_venv(ctx, kwds)
        run_script = "%s $COMMON_STARTUP_ARGS" % shlex_quote(os.path.join(self.galaxy_root, "run.sh"))
        if daemon:
            run_script += " --daemon"
            self.env["GALAXY_RUN_ALL"] = "1"
        else:
            run_script += " --server-name %s" % shlex_quote(self.server_name)
        server_ini = os.path.join(self.config_directory, "galaxy.ini")
        self.env["GALAXY_CONFIG_FILE"] = server_ini
        if parse_version(kwds.get('galaxy_python_version') or DEFAULT_PYTHON_VERSION) >= parse_version('3'):
            # We need to start under gunicorn
            self.env['APP_WEBSERVER'] = 'gunicorn'
            self.env['GUNICORN_CMD_ARGS'] = "--timeout={timeout} --capture-output --bind={host}:{port} --name={server_name}".format(
                timeout=DEFAULT_TOOL_INSTALL_TIMEOUT,
                host=kwds.get('host', '127.0.0.1'),
                port=kwds['port'],
                server_name=self.server_name,
            )
        cd_to_galaxy_command = ['cd', self.galaxy_root]
        return shell_join(
            cd_to_galaxy_command,
            setup_venv_command,
            setup_common_startup_args(),
            run_script,
        )

    @property
    def log_file(self):
        """Log file used when planemo serves this Galaxy instance."""
        file_name = "%s.log" % self.server_name
        return os.path.join(self.galaxy_root, file_name)

    @property
    def pid_file(self):
        pid_file_name = "%s.pid" % self.server_name
        return os.path.join(self.galaxy_root, pid_file_name)

    @property
    def log_contents(self):
        if not os.path.exists(self.log_file):
            return ""
        with open(self.log_file, "r") as f:
            return f.read()

    def cleanup(self):
        shutil.rmtree(self.config_directory, CLEANUP_IGNORE_ERRORS)

    @property
    def default_use_path_paste(self):
        # If Planemo started a local, native Galaxy instance assume files URLs can be
        # pasted.
        return self.user_is_admin


def _database_connection(database_location, **kwds):
    default_connection = DATABASE_LOCATION_TEMPLATE % database_location
    database_connection = kwds.get("database_connection") or default_connection
    return database_connection


def _find_galaxy_root(ctx, **kwds):
    root_prop = "galaxy_root"
    cwl = kwds.get("cwl", False)
    if cwl:
        root_prop = "cwl_galaxy_root"
    galaxy_root = kwds.get(root_prop, None)
    if galaxy_root:
        return galaxy_root
    else:
        par_dir = os.getcwd()
        while True:
            run = os.path.join(par_dir, "run.sh")
            config = os.path.join(par_dir, "config")
            if os.path.isfile(run) and os.path.isdir(config):
                return par_dir
            new_par_dir = os.path.dirname(par_dir)
            if new_par_dir == par_dir:
                break
            par_dir = new_par_dir
    return None


def _find_test_data(runnables, **kwds):
    test_data_search_path = "."
    runnables = [r for r in runnables if r.has_tools]
    if len(runnables) > 0:
        test_data_search_path = runnables[0].test_data_search_path

    # Find test data directory associated with path.
    test_data = kwds.get("test_data", None)
    if test_data:
        return os.path.abspath(test_data)
    else:
        test_data = _search_tool_path_for(test_data_search_path, "test-data")
        if test_data:
            return test_data
    warn(NO_TEST_DATA_MESSAGE)
    return None


def _find_tool_data_table(runnables, test_data_dir, **kwds):
    tool_data_search_path = "."
    runnables = [r for r in runnables if r.has_tools]
    if len(runnables) > 0:
        tool_data_search_path = runnables[0].tool_data_search_path

    tool_data_table = kwds.get("tool_data_table", None)
    if tool_data_table:
        return os.path.abspath(tool_data_table)
    else:
        extra_paths = [test_data_dir] if test_data_dir else []
        return _search_tool_path_for(
            tool_data_search_path,
            "tool_data_table_conf.xml.test",
            extra_paths,
        ) or _search_tool_path_for(  # if all else fails just use sample
            tool_data_search_path,
            "tool_data_table_conf.xml.sample"
        )


def _search_tool_path_for(path, target, extra_paths=None):
    """Check for presence of a target in different artifact directories."""
    if extra_paths is None:
        extra_paths = []
    if not os.path.isdir(path):
        tool_dir = os.path.dirname(path)
    else:
        tool_dir = path
    possible_dirs = [tool_dir, "."] + extra_paths
    for possible_dir in possible_dirs:
        possible_path = os.path.join(possible_dir, target)
        if os.path.exists(possible_path):
            return os.path.abspath(possible_path)
    return None


def _configure_sheds_config_file(ctx, config_directory, **kwds):
    if "shed_target" not in kwds:
        kwds = kwds.copy()
        kwds["shed_target"] = "toolshed"
    shed_target_url = tool_shed_url(ctx, **kwds)
    contents = _sub(TOOL_SHEDS_CONF, {"shed_target_url": shed_target_url})
    tool_sheds_conf = os.path.join(config_directory, "tool_sheds_conf.xml")
    write_file(tool_sheds_conf, contents)
    return tool_sheds_conf


def _tool_conf_entry_for(tool_paths):
    tool_definitions = ""
    for tool_path in tool_paths:
        if os.path.isdir(tool_path):
            tool_definitions += '''<tool_dir dir="%s" />''' % tool_path
        else:
            tool_definitions += '''<tool file="%s" />''' % tool_path
    return tool_definitions


def _install_galaxy(ctx, galaxy_root, env, kwds):
    if not kwds.get("no_cache_galaxy", False):
        _install_galaxy_via_git(ctx, galaxy_root, env, kwds)
    else:
        _install_galaxy_via_download(ctx, galaxy_root, env, kwds)


def _install_galaxy_via_download(ctx, galaxy_root, env, kwds):
    branch = _galaxy_branch(kwds)
    untar_to("https://codeload.github.com/galaxyproject/galaxy/tar.gz/" + branch, tar_args=['-xvzf', '-', 'galaxy-' + branch], dest_dir=galaxy_root)
    _install_with_command(ctx, galaxy_root, env, kwds)


def _install_galaxy_via_git(ctx, galaxy_root, env, kwds):
    gx_repo = _ensure_galaxy_repository_available(ctx, kwds)
    branch = _galaxy_branch(kwds)
    command = git.command_clone(ctx, gx_repo, galaxy_root, branch=branch)
    exit_code = shell(command, env=env)
    if exit_code != 0:
        raise Exception("Failed to glone Galaxy via git")
    _install_with_command(ctx, galaxy_root, env, kwds)


def _build_eggs_cache(ctx, env, kwds):
    if kwds.get("no_cache_galaxy", False):
        return None
    workspace = ctx.workspace
    eggs_path = os.path.join(workspace, "gx_eggs")
    if not os.path.exists(eggs_path):
        os.makedirs(eggs_path)
    env["GALAXY_EGGS_PATH"] = eggs_path


def _galaxy_branch(kwds):
    branch = kwds.get("galaxy_branch", None)
    if branch is None:
        cwl = kwds.get("cwl", False)
        branch = "cwl-1.0" if cwl else None
    if branch is None:
        branch = DEFAULT_GALAXY_BRANCH

    return branch


def _galaxy_source(kwds):
    source = kwds.get("galaxy_source", None)
    if source is None:
        cwl = kwds.get("cwl", False)
        source = CWL_GALAXY_SOURCE if cwl else None
    if source is None:
        source = DEFAULT_GALAXY_SOURCE

    return source


def _install_with_command(ctx, galaxy_root, env, kwds):
    setup_venv_command = setup_venv(ctx, kwds)
    env['__PYVENV_LAUNCHER__'] = ''
    install_cmd = shell_join(
        setup_venv_command,
        setup_common_startup_args(),
        COMMAND_STARTUP_COMMAND,
    )
    exit_code = shell(install_cmd, cwd=galaxy_root, env=env)
    if exit_code != 0:
        raise Exception("Failed to install Galaxy via command [%s]" % install_cmd)
    if not os.path.exists(galaxy_root):
        raise Exception("Failed to create Galaxy directory [%s]" % galaxy_root)
    if not os.path.exists(os.path.join(galaxy_root, "lib")):
        raise Exception("Failed to create Galaxy directory [%s], lib missing" % galaxy_root)


def _ensure_galaxy_repository_available(ctx, kwds):
    workspace = ctx.workspace
    cwl = kwds.get("cwl", False)
    galaxy_source = kwds.get('galaxy_source')
    if galaxy_source and galaxy_source != DEFAULT_GALAXY_SOURCE:
        sanitized_repo_name = "".join(c if c.isalnum() else '_' for c in kwds['galaxy_source']).rstrip()[:255]
        gx_repo = os.path.join(workspace, "gx_repo_%s" % sanitized_repo_name)
    else:
        gx_repo = os.path.join(workspace, "gx_repo")
    if cwl:
        gx_repo += "_cwl"
    if os.path.exists(gx_repo):
        # Convert the git repository from bare to mirror, if needed
        shell(['git', '--git-dir', gx_repo, 'config', 'remote.origin.fetch', '+refs/*:refs/*'])
        shell(['git', '--git-dir', gx_repo, 'config', 'remote.origin.mirror', 'true'])
        # Attempt remote update - but don't fail if not interweb, etc...
        shell("git --git-dir %s remote update >/dev/null 2>&1" % gx_repo)
    else:
        remote_repo = _galaxy_source(kwds)
        command = git.command_clone(ctx, remote_repo, gx_repo, mirror=True)
        shell(command)
    return gx_repo


def _build_env_for_galaxy(properties, template_args):
    env = {}
    for key, value in iteritems(properties):
        if value is not None:  # Do not override None with empty string
            var = "GALAXY_CONFIG_OVERRIDE_%s" % key.upper()
            value = _sub(value, template_args)
            env[var] = value
    return env


def _build_test_env(properties, env):
    # Keeping these environment variables around for a little while but
    # many are probably not needed as of the following commit.
    # https://bitbucket.org/galaxy/galaxy-central/commits/d7dd1f9
    test_property_variants = {
        'GALAXY_TEST_JOB_CONFIG_FILE': 'job_config_file',
        'GALAXY_TEST_MIGRATED_TOOL_CONF': 'migrated_tools_config',
        'GALAXY_TEST_TOOL_CONF': 'tool_config_file',
        'GALAXY_TEST_FILE_DIR': 'test_data_dir',
        'GALAXY_TOOL_DEPENDENCY_DIR': 'tool_dependency_dir',
        # Next line would be required for tool shed tests.
        # 'GALAXY_TEST_TOOL_DEPENDENCY_DIR': 'tool_dependency_dir',
    }
    for test_key, gx_key in test_property_variants.items():
        value = properties.get(gx_key, None)
        if value is not None:
            env[test_key] = value


def _handle_job_config_file(config_directory, server_name, kwds):
    job_config_file = kwds.get("job_config_file", None)
    if not job_config_file:
        template_str = JOB_CONFIG_LOCAL
        job_config_file = os.path.join(
            config_directory,
            "job_conf.xml",
        )
        docker_enable = str(kwds.get("docker", False))
        docker_host = str(kwds.get("docker_host", docker_util.DEFAULT_HOST))
        docker_host_param = ""
        if docker_host:
            docker_host_param = """<param id="docker_host">%s</param>""" % docker_host

        conf_contents = Template(template_str).safe_substitute({
            "server_name": server_name,
            "docker_enable": docker_enable,
            "require_container": "false",
            "docker_sudo": str(kwds.get("docker_sudo", False)),
            "docker_sudo_cmd": str(kwds.get("docker_sudo_cmd", docker_util.DEFAULT_SUDO_COMMAND)),
            "docker_cmd": str(kwds.get("docker_cmd", docker_util.DEFAULT_DOCKER_COMMAND)),
            "docker_host": docker_host_param,
        })
        write_file(job_config_file, conf_contents)
    kwds["job_config_file"] = job_config_file


def _write_tool_conf(ctx, tool_paths, tool_conf_path):
    tool_definition = _tool_conf_entry_for(tool_paths)
    tool_conf_template_kwds = dict(tool_definition=tool_definition)
    tool_conf_contents = _sub(TOOL_CONF_TEMPLATE, tool_conf_template_kwds)
    write_file(tool_conf_path, tool_conf_contents)
    ctx.vlog(
        "Writing tool_conf to path %s with contents [%s]",
        tool_conf_path,
        tool_conf_contents,
    )


def _handle_container_resolution(ctx, kwds, galaxy_properties):
    if kwds.get("mulled_containers", False):
        galaxy_properties["enable_beta_mulled_containers"] = "True"
        involucro_context = build_involucro_context(ctx, **kwds)
        galaxy_properties["involucro_auto_init"] = "False"  # Use planemo's
        galaxy_properties["involucro_path"] = involucro_context.involucro_bin


def _handle_job_metrics(config_directory, kwds):
    metrics_conf = os.path.join(config_directory, "job_metrics_conf.xml")
    with open(metrics_conf, "w") as fh:
        fh.write(EMPTY_JOB_METRICS_TEMPLATE)
    kwds["job_metrics_config_file"] = metrics_conf


def _handle_refgenie_config(config_directory, kwds):
    refgenie_dir = os.path.join(config_directory, 'refgenie')
    _ensure_directory(refgenie_dir)
    refgenie_config = os.path.join(refgenie_dir, "genome_config.yaml")
    with open(refgenie_config, "w") as fh:
        fh.write(REFGENIE_CONFIG_TEMPLATE % (refgenie_dir))
    kwds["refgenie_config_file"] = refgenie_config


def _handle_kwd_overrides(properties, kwds):
    kwds_gx_properties = [
        'job_config_file',
        'job_metrics_config_file',
        'dependency_resolvers_config_file',
    ]
    for prop in kwds_gx_properties:
        val = kwds.get(prop, None)
        if val:
            properties[prop] = val


def _sub(template, args):
    if template is None:
        return ''
    return Template(template).safe_substitute(args)


def _ensure_directory(path):
    if path is not None and not os.path.exists(path):
        os.makedirs(path)


__all__ = (
    "DATABASE_LOCATION_TEMPLATE",
    "galaxy_config",
)
