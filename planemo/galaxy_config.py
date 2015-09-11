from __future__ import absolute_import
from __future__ import print_function

import contextlib
import os
import random
import shutil
import time
from six.moves.urllib.request import urlopen
from six import iteritems
from string import Template
from tempfile import mkdtemp
from six.moves.urllib.request import urlretrieve

import click

from planemo import galaxy_run
from planemo.io import warn
from planemo.io import shell
from planemo.io import write_file
from planemo.io import kill_pid_file
from planemo import git
from planemo.shed import tool_shed_url
from planemo.bioblend import (
    galaxy,
    ensure_module,
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
"""

TOOL_CONF_TEMPLATE = """<toolbox>
  <tool file="data_source/upload.xml" />
  ${tool_definition}
</toolbox>
"""

SHED_TOOL_CONF_TEMPLATE = """<?xml version="1.0"?>
<toolbox tool_path="${shed_tools_path}">
</toolbox>
"""


EMPTY_JOB_METRICS_TEMPLATE = """<?xml version="1.0"?>
<job_metrics>
</job_metrics>
"""


BREW_DEPENDENCY_RESOLUTION_CONF = """<dependency_resolvers>
  <homebrew />
  <!--
  <homebrew versionless="true" />
  -->
</dependency_resolvers>
"""

SHED_DEPENDENCY_RESOLUTION_CONF = """<dependency_resolvers>
  <tool_shed_tap />
</dependency_resolvers>
"""

TOOL_SHEDS_CONF = """<tool_sheds>
  <tool_shed name="Target Shed" url="${shed_target_url}" />
</tool_sheds>
"""

# Provide some shortcuts for simple/common dependency resolutions strategies.
STOCK_DEPENDENCY_RESOLUTION_STRATEGIES = {
    "brew_dependency_resolution": BREW_DEPENDENCY_RESOLUTION_CONF,
    "shed_dependency_resolution": SHED_DEPENDENCY_RESOLUTION_CONF,
}

EMPTY_TOOL_CONF_TEMPLATE = """<toolbox></toolbox>"""

DOWNLOADS_URL = ("https://raw.githubusercontent.com/"
                 "jmchilton/galaxy-downloads/master/")
DOWNLOADABLE_MIGRATION_VERSIONS = [127, 120, 117]
LATEST_URL = DOWNLOADS_URL + "latest.sqlite"

FAILED_TO_FIND_GALAXY_EXCEPTION = (
    "Failed to find Galaxy root directory - please explicitly specify one "
    "with --galaxy_root."
)


@contextlib.contextmanager
def galaxy_config(ctx, tool_paths, for_tests=False, **kwds):
    test_data_dir = _find_test_data(tool_paths, **kwds)
    tool_data_table = _find_tool_data_table(
        tool_paths,
        test_data_dir=test_data_dir,
        **kwds
    )
    galaxy_root = _check_galaxy(ctx, **kwds)
    install_galaxy = galaxy_root is None
    config_directory = kwds.get("config_directory", None)

    def config_join(*args):
        return os.path.join(config_directory, *args)

    created_config_directory = False
    if not config_directory:
        created_config_directory = True
        config_directory = mkdtemp()
    try:
        latest_galaxy = False
        if install_galaxy:
            _install_galaxy(ctx, config_directory, kwds)
            latest_galaxy = True
            galaxy_root = config_join("galaxy-dev")

        _handle_dependency_resolution(config_directory, kwds)
        _handle_job_metrics(config_directory, kwds)
        tool_definition = _tool_conf_entry_for(tool_paths)
        empty_tool_conf = config_join("empty_tool_conf.xml")
        shed_tool_conf = _shed_tool_conf(install_galaxy, config_directory)
        tool_conf = config_join("tool_conf.xml")
        database_location = config_join("galaxy.sqlite")
        shed_tools_path = config_join("shed_tools")
        sheds_config_path = _configure_sheds_config_file(
            config_directory, **kwds
        )
        preseeded_database = True
        master_api_key = kwds.get("master_api_key", "test_key")
        dependency_dir = os.path.join(config_directory, "deps")

        try:
            _download_database_template(
                galaxy_root,
                database_location,
                latest=latest_galaxy
            )
        except Exception as e:
            print(e)
            # No network access - just roll forward from null.
            preseeded_database = False

        os.makedirs(shed_tools_path)
        server_name = "planemo%d" % random.randint(0, 100000)
        port = kwds.get("port", 9090)
        template_args = dict(
            port=port,
            host=kwds.get("host", "127.0.0.1"),
            server_name=server_name,
            temp_directory=config_directory,
            shed_tools_path=shed_tools_path,
            database_location=database_location,
            tool_definition=tool_definition,
            tool_conf=tool_conf,
            debug=kwds.get("debug", "true"),
            master_api_key=master_api_key,
            id_secret=kwds.get("id_secret", "test_secret"),
            log_level=kwds.get("log_level", "DEBUG"),
        )
        tool_config_file = "%s,%s" % (tool_conf, shed_tool_conf)
        properties = dict(
            tool_dependency_dir=dependency_dir,
            file_path="${temp_directory}/files",
            new_file_path="${temp_directory}/tmp",
            tool_config_file=tool_config_file,
            tool_sheds_config_file=sheds_config_path,
            check_migrate_tools="False",
            manage_dependency_relationships="False",
            job_working_directory="${temp_directory}/job_working_directory",
            template_cache_path="${temp_directory}/compiled_templates",
            citation_cache_type="file",
            citation_cache_data_dir="${temp_directory}/citations/data",
            citation_cache_lock_dir="${temp_directory}/citations/lock",
            collect_outputs_from="job_working_directory",
            database_auto_migrate="True",
            cleanup_job="never",
            master_api_key="${master_api_key}",
            id_secret="${id_secret}",
            log_level="${log_level}",
            debug="${debug}",
            watch_tools="auto",
            tool_data_table_config_path=tool_data_table,
            integrated_tool_panel_config=("${temp_directory}/"
                                          "integrated_tool_panel_conf.xml"),
            # Use in-memory database for kombu to avoid database contention
            # during tests.
            amqp_internal_connection="sqlalchemy+sqlite://",
            migrated_tools_config=empty_tool_conf,
            test_data_dir=test_data_dir,  # TODO: make gx respect this
        )
        if not for_tests:
            properties["database_connection"] = \
                "sqlite:///${database_location}?isolation_level=IMMEDIATE"

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
        if install_galaxy:
            _build_eggs_cache(ctx, env, kwds)
        _build_test_env(properties, env)
        env['GALAXY_TEST_SHED_TOOL_CONF'] = shed_tool_conf

        # No need to download twice - would GALAXY_TEST_DATABASE_CONNECTION
        # work?
        if preseeded_database:
            env["GALAXY_TEST_DB_TEMPLATE"] = os.path.abspath(database_location)
        env["GALAXY_TEST_UPLOAD_ASYNC"] = "false"
        env["GALAXY_DEVELOPMENT_ENVIRONMENT"] = "1"
        web_config = _sub(WEB_SERVER_CONFIG_TEMPLATE, template_args)
        write_file(config_join("galaxy.ini"), web_config)
        tool_conf_contents = _sub(TOOL_CONF_TEMPLATE, template_args)
        write_file(tool_conf, tool_conf_contents)
        write_file(empty_tool_conf, EMPTY_TOOL_CONF_TEMPLATE)

        shed_tool_conf_contents = _sub(SHED_TOOL_CONF_TEMPLATE, template_args)
        write_file(shed_tool_conf, shed_tool_conf_contents)

        yield GalaxyConfig(
            galaxy_root,
            config_directory,
            env,
            test_data_dir,
            port,
            server_name,
            master_api_key,
        )
    finally:
        cleanup = not kwds.get("no_cleanup", False)
        if created_config_directory and cleanup:
            shutil.rmtree(config_directory)


class GalaxyConfig(object):

    def __init__(
        self,
        galaxy_root,
        config_directory,
        env,
        test_data_dir,
        port,
        server_name,
        master_api_key,
    ):
        self.galaxy_root = galaxy_root
        self.config_directory = config_directory
        self.env = env
        self.test_data_dir = test_data_dir
        # Runtime server configuration stuff not used if testing...
        # better design might be GalaxyRootConfig and GalaxyServerConfig
        # as two separate objects.
        self.port = port
        self.server_name = server_name
        self.master_api_key = master_api_key

    def kill(self):
        kill_pid_file(self.pid_file)

    @property
    def pid_file(self):
        return os.path.join(self.galaxy_root, "%s.pid" % self.server_name)

    @property
    def gi(self):
        ensure_module(galaxy)
        return galaxy.GalaxyInstance(
            url="http://localhost:%d" % self.port,
            key=self.master_api_key
        )

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
                return False
            if status == "Installed":
                return True
            raise Exception("Error installing repo status is %s" % status)

        def not_ready():
            repos = self.tool_shed_client.get_repositories()
            return not all(map(status_ready, repos))

        self._wait_for(not_ready)

    # Taken from Galaxy's twilltestcase.
    def _wait_for(self, func, **kwd):
        sleep_amount = 0.2
        slept = 0
        walltime_exceeded = 1086400
        while slept <= walltime_exceeded:
            result = func()
            if result:
                time.sleep(sleep_amount)
                slept += sleep_amount
                sleep_amount *= 1.25
                if slept + sleep_amount > walltime_exceeded:
                    sleep_amount = walltime_exceeded - slept
            else:
                break
        assert slept < walltime_exceeded, "Action taking too long."

    def cleanup(self):
        shutil.rmtree(self.config_directory)


def _download_database_template(galaxy_root, database_location, latest=False):
    if latest:
        template_url = DOWNLOADS_URL + urlopen(LATEST_URL).read()
        urlretrieve(template_url, database_location)
        return True

    newest_migration = _newest_migration_version(galaxy_root)
    download_migration = None
    for migration in DOWNLOADABLE_MIGRATION_VERSIONS:
        if newest_migration > migration:
            download_migration = migration
            break

    if download_migration:
        download_name = "db_gx_rev_0%d.sqlite" % download_migration
        download_url = DOWNLOADS_URL + download_name
        urlretrieve(download_url, database_location)
        return True
    else:
        return False


def _newest_migration_version(galaxy_root):
    versions = os.path.join(galaxy_root, "lib/galaxy/model/migrate/versions")
    version = max(map(_file_name_to_migration_version, os.listdir(versions)))
    return version


def _file_name_to_migration_version(name):
    try:
        return int(name[0:4])
    except ValueError:
        return -1


def _check_galaxy(ctx, **kwds):
    """ Find Galaxy root, return None to indicate it should be
    installed automatically.
    """
    install_galaxy = kwds.get("install_galaxy", None)
    galaxy_root = None
    if not install_galaxy:
        galaxy_root = _find_galaxy_root(ctx, **kwds)
    return galaxy_root


def _find_galaxy_root(ctx, **kwds):
    galaxy_root = kwds.get("galaxy_root", None)
    if galaxy_root:
        return galaxy_root
    elif ctx.global_config.get("galaxy_root", None):
        return ctx.global_config["galaxy_root"]
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


def _find_test_data(tool_paths, **kwds):
    path = "."
    if len(tool_paths) > 0:
        path = tool_paths[0]

    # Find test data directory associated with path.
    test_data = kwds.get("test_data", None)
    if test_data:
        return os.path.abspath(test_data)
    else:
        test_data = _search_tool_path_for(path, "test-data")
        if test_data:
            return test_data
    warn(NO_TEST_DATA_MESSAGE)
    return None


def _find_tool_data_table(tool_paths, test_data_dir, **kwds):
    path = "."
    if len(tool_paths) > 0:
        path = tool_paths[0]

    tool_data_table = kwds.get("tool_data_table", None)
    if tool_data_table:
        return os.path.abspath(tool_data_table)
    else:
        extra_paths = [test_data_dir] if test_data_dir else []
        return _search_tool_path_for(
            path,
            "tool_data_table_conf.xml.test",
            extra_paths,
        ) or _search_tool_path_for(  # if all else fails just use sample
            path,
            "tool_data_table_conf.xml.sample"
        )


def _search_tool_path_for(path, target, extra_paths=[]):
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


def _configure_sheds_config_file(config_directory, **kwds):
    if "shed_target" not in kwds:
        kwds = kwds.copy()
        kwds["shed_target"] = "toolshed"
    shed_target_url = tool_shed_url(kwds)
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


def _shed_tool_conf(install_galaxy, config_directory):
    # TODO: There is probably a reason this is split up like this but I have
    # no clue why I did it and not documented on the commit message.
    if install_galaxy:
        config_dir = os.path.join(config_directory, "galaxy-dev", "config")
    else:
        config_dir = config_directory
    return os.path.join(config_dir, "shed_tool_conf.xml")


def _install_galaxy(ctx, config_directory, kwds):
    if not kwds.get("no_cache_galaxy", False):
        _install_galaxy_via_git(ctx, config_directory, kwds)
    else:
        _install_galaxy_via_download(config_directory, kwds)


def _install_galaxy_via_download(config_directory, kwds):
    command = galaxy_run.DOWNLOAD_GALAXY + "; tar -zxvf dev | tail"
    _install_with_command(config_directory, command)


def _install_galaxy_via_git(ctx, config_directory, kwds):
    _ensure_galaxy_repository_available(ctx)
    workspace = ctx.workspace
    gx_repo = os.path.join(workspace, "gx_repo")
    command = git.command_clone(ctx, gx_repo, "galaxy-dev")
    _install_with_command(config_directory, command)


def _build_eggs_cache(ctx, env, kwds):
    if kwds.get("no_cache_galaxy", False):
        return None
    workspace = ctx.workspace
    eggs_path = os.path.join(workspace, "gx_eggs")
    if not os.path.exists(eggs_path):
        os.makedirs(eggs_path)
    env["GALAXY_EGGS_PATH"] = eggs_path


def _install_with_command(config_directory, command):
    install_cmds = [
        "cd %s" % config_directory,
        command,
        "cd galaxy-dev",
        "type virtualenv >/dev/null 2>&1 && virtualenv .venv",
        galaxy_run.ACTIVATE_COMMAND,
    ]
    shell(";".join(install_cmds))


def _ensure_galaxy_repository_available(ctx):
    workspace = ctx.workspace
    gx_repo = os.path.join(workspace, "gx_repo")
    if os.path.exists(gx_repo):
        # Attempt fetch - but don't fail if not interweb, etc...
        shell("git --git-dir %s fetch >/dev/null 2>&1" % gx_repo)
    else:
        remote_repo = "https://github.com/galaxyproject/galaxy"
        command = git.command_clone(ctx, remote_repo, gx_repo, bare=True)
        shell(command)


def _build_env_for_galaxy(properties, template_args):
    env = {}
    for key, value in iteritems(properties):
        var = "GALAXY_CONFIG_OVERRIDE_%s" % key.upper()
        value = _sub(value, template_args)
        env[var] = value
    return env


def _build_test_env(properties, env):
    # Keeping these environment variables around for a little while but they
    # many are probably not needed as of the following commit.
    # https://bitbucket.org/galaxy/galaxy-central/commits/d7dd1f9
    test_property_variants = {
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


def _handle_dependency_resolution(config_directory, kwds):
    resolutions_strategies = [
        "brew_dependency_resolution",
        "dependency_resolvers_config_file",
        "shed_dependency_resolution",
    ]

    selected_strategies = 0
    for key in resolutions_strategies:
        if kwds.get(key):
            selected_strategies += 1

    if selected_strategies > 1:
        message = "At most one option from [%s] may be specified"
        raise click.UsageError(message % resolutions_strategies)

    for key in STOCK_DEPENDENCY_RESOLUTION_STRATEGIES:
        if kwds.get(key):
            resolvers_conf = os.path.join(
                config_directory,
                "resolvers_conf.xml"
            )
            conf_contents = STOCK_DEPENDENCY_RESOLUTION_STRATEGIES[key]
            open(resolvers_conf, "w").write(conf_contents)
            kwds["dependency_resolvers_config_file"] = resolvers_conf


def _handle_job_metrics(config_directory, kwds):
    metrics_conf = os.path.join(config_directory, "job_metrics_conf.xml")
    open(metrics_conf, "w").write(EMPTY_JOB_METRICS_TEMPLATE)
    kwds["job_metrics_config_file"] = metrics_conf


def _handle_kwd_overrides(properties, kwds):
    kwds_gx_properties = [
        'job_config_file',
        'job_metrics_config_file',
        'dependency_resolvers_config_file',
        'tool_dependency_dir',
    ]
    for prop in kwds_gx_properties:
        val = kwds.get(prop, None)
        if val:
            properties[prop] = val


def _sub(template, args):
    if template is None:
        return ''
    return Template(template).safe_substitute(args)

__all__ = ["galaxy_config"]
