from __future__ import absolute_import
from collections import namedtuple
import contextlib
import os
import shutil
import urllib
from string import Template
from tempfile import mkdtemp
from six.moves.urllib.request import urlretrieve

import click

from planemo import galaxy_run
from planemo.io import warn

from galaxy.tools.deps import commands

NO_TEST_DATA_MESSAGE = (
    "planemo couldn't find a target test-data directory, you should likely "
    "create a test-data directory of pass an explicit path using --test-data."
)

WEB_SERVER_CONFIG_TEMPLATE = """
[server:main]
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

GalaxyConfig = namedtuple(
    'GalaxyConfig',
    ['galaxy_root', 'config_directory', 'env', 'test_data_dir']
)


@contextlib.contextmanager
def galaxy_config(ctx, tool_path, for_tests=False, **kwds):
    test_data_dir = __find_test_data(tool_path, **kwds)
    tool_data_table = __find_tool_data_table(
        tool_path,
        test_data_dir=test_data_dir,
        **kwds
    )
    if kwds.get("install_galaxy", None):
        galaxy_root = None
    else:
        galaxy_root = __find_galaxy_root(ctx, **kwds)

    config_directory = kwds.get("config_directory", None)

    def config_join(*args):
        return os.path.join(config_directory, *args)

    created_config_directory = False
    if not config_directory:
        created_config_directory = True
        config_directory = mkdtemp()
    try:
        latest_galaxy = False
        if __install_galaxy_if_needed(config_directory, kwds):
            latest_galaxy = True
            galaxy_root = config_join("galaxy-central-master")

        __handle_dependency_resolution(config_directory, kwds)
        __handle_job_metrics(config_directory, kwds)
        tool_definition = __tool_conf_entry_for(tool_path)
        empty_tool_conf = config_join("empty_tool_conf.xml")
        tool_conf = config_join("tool_conf.xml")
        database_location = config_join("galaxy.sqlite")
        preseeded_database = True

        try:
            _download_database_template(
                galaxy_root,
                database_location,
                latest=latest_galaxy
            )
        except Exception as e:
            print e
            # No network access - just roll forward from null.
            preseeded_database = False

        template_args = dict(
            port=kwds.get("port", 9090),
            host="127.0.0.1",
            temp_directory=config_directory,
            database_location=database_location,
            tool_definition=tool_definition % tool_path,
            tool_conf=tool_conf,
            debug=kwds.get("debug", "true"),
            master_api_key=kwds.get("master_api_key", "test_key"),
            id_secret=kwds.get("id_secret", "test_secret"),
            log_level=kwds.get("log_level", "DEBUG"),
        )
        properties = dict(
            file_path="${temp_directory}files",
            new_file_path="${temp_directory}/tmp",
            tool_config_file=tool_conf,
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

        __handle_kwd_overrides(properties, kwds)

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

        env = __build_env_for_galaxy(properties, template_args)
        __build_test_env(properties, env)

        # No need to download twice - would GALAXY_TEST_DATABASE_CONNECTION
        # work?
        if preseeded_database:
            env["GALAXY_TEST_DB_TEMPLATE"] = os.path.abspath(database_location)
        env["GALAXY_TEST_UPLOAD_ASYNC"] = "false"
        env["GALAXY_DEVELOPMENT_ENVIRONMENT"] = "1"
        web_config = __sub(WEB_SERVER_CONFIG_TEMPLATE, template_args)
        open(config_join("galaxy.ini"), "w").write(web_config)
        tool_conf_contents = __sub(TOOL_CONF_TEMPLATE, template_args)
        open(tool_conf, "w").write(tool_conf_contents)
        open(empty_tool_conf, "w").write(EMPTY_TOOL_CONF_TEMPLATE)

        yield GalaxyConfig(galaxy_root, config_directory, env, test_data_dir)
    finally:
        if created_config_directory:
            shutil.rmtree(config_directory)


def _download_database_template(galaxy_root, database_location, latest=False):
    if latest:
        template_url = DOWNLOADS_URL + urllib.urlopen(LATEST_URL).read()
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


def __find_galaxy_root(ctx, **kwds):
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
    raise Exception(FAILED_TO_FIND_GALAXY_EXCEPTION)


def __find_test_data(path, **kwds):
    # Find test data directory associated with path.
    test_data = kwds.get("test_data", None)
    if test_data:
        return os.path.abspath(test_data)
    else:
        test_data = __search_tool_path_for(path, "test-data")
        if test_data:
            return test_data
    warn(NO_TEST_DATA_MESSAGE)
    return None


def __find_tool_data_table(path, test_data_dir, **kwds):
    tool_data_table = kwds.get("tool_data_table", None)
    if tool_data_table:
        return os.path.abspath(tool_data_table)
    else:
        extra_paths = [test_data_dir] if test_data_dir else []
        return __search_tool_path_for(
            path,
            "tool_data_table_conf.xml.test",
            extra_paths,
        ) or __search_tool_path_for(  # if all else fails just use sample
            path,
            "tool_data_table_conf.xml.sample"
        )


def __search_tool_path_for(path, target, extra_paths=[]):
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


def __tool_conf_entry_for(tool_path):
    if os.path.isdir(tool_path):
        tool_definition = '''<tool_dir dir="%s" />'''
    else:
        tool_definition = '''<tool file="%s" />'''
    return tool_definition


def __install_galaxy_if_needed(config_directory, kwds):
    installed = False
    if kwds.get("install_galaxy", None):
        install_cmds = [
            "cd %s" % config_directory,
            galaxy_run.DOWNLOAD_GALAXY,
            "tar -zxvf master | tail",
            "cd galaxy-central-master",
            "type virtualenv >/dev/null 2>&1 && virtualenv .venv",
            galaxy_run.ACTIVATE_COMMAND,
        ]
        commands.shell(";".join(install_cmds))
        installed = True
    return installed


def __build_env_for_galaxy(properties, template_args):
    env = {}
    for key, value in properties.iteritems():
        var = "GALAXY_CONFIG_OVERRIDE_%s" % key.upper()
        value = __sub(value, template_args)
        env[var] = value
    return env


def __build_test_env(properties, env):
    # Keeping these environment variables around for a little while but they
    # many are probably not needed as of the following commit.
    # https://bitbucket.org/galaxy/galaxy-central/commits/d7dd1f9
    test_property_variants = {
        'GALAXY_TEST_MIGRATED_TOOL_CONF': 'migrated_tools_config',
        'GALAXY_TEST_SHED_TOOL_CONF': 'migrated_tools_config',  # Hack
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


def __handle_dependency_resolution(config_directory, kwds):
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
            dependency_dir = os.path.join("config_directory", "deps")
            kwds["tool_dependency_dir"] = dependency_dir
            kwds["dependency_resolvers_config_file"] = resolvers_conf


def __handle_job_metrics(config_directory, kwds):
    metrics_conf = os.path.join(config_directory, "job_metrics_conf.xml")
    open(metrics_conf, "w").write(EMPTY_JOB_METRICS_TEMPLATE)
    kwds["job_metrics_config_file"] = metrics_conf


def __handle_kwd_overrides(properties, kwds):
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


def __sub(template, args):
    if template is None:
        return ''
    return Template(template).safe_substitute(args)

__all__ = ["galaxy_config"]
