"""Unit tests for ``planemo.galaxy.config``."""
import contextlib
import os

import yaml

from planemo.galaxy.config import (
    galaxy_config,
    get_refgenie_config,
)
from .test_utils import (
    create_test_context,
    skip_if_environ,
    TempDirectoryContext,
)


@skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
def test_defaults():
    """Test by default Galaxy files are stored in temp ``config_directory``."""
    with _test_galaxy_config() as config:
        config_directory = config.config_directory
        _assert_property_is(config, "file_path", os.path.join(config_directory, "files"))
        conn = "sqlite:///%s/galaxy.sqlite?isolation_level=IMMEDIATE" % config_directory
        _assert_property_is(config, "database_connection", conn)


@skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
def test_database_connection_override_path():
    """Test by default Galaxy files are stored in temp ``config_directory``."""
    conn = "postgresql://username:password@localhost/mydatabase"
    with _test_galaxy_config(database_connection=conn) as config:
        _assert_property_is(config, "database_connection", conn)


@skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
def test_override_files_path():
    """Test Galaxy file path overrideable with --file_path."""
    with TempDirectoryContext() as tdc:
        with _test_galaxy_config(file_path=tdc.temp_directory) as config:
            _assert_property_is(config, "file_path", tdc.temp_directory)


def test_refgenie_config_version():
    with TempDirectoryContext() as tdc:
        galaxy_lib_path = os.path.join(tdc.temp_directory, "lib", "galaxy")
        os.makedirs(galaxy_lib_path)
        version_path = os.path.join(galaxy_lib_path, "version.py")
        with open(version_path, "w") as version_fh:
            version_fh.write('VERSION_MAJOR = "21.05"')
        refgenie_config = get_refgenie_config(galaxy_root=tdc.temp_directory, refgenie_dir="/")
    assert yaml.load(refgenie_config, Loader=yaml.SafeLoader)["config_version"] == 0.3


def _assert_property_is(config, prop, value):
    env_var = "GALAXY_CONFIG_OVERRIDE_%s" % prop.upper()
    assert config.env[env_var] == value


@contextlib.contextmanager
def _test_galaxy_config(tool_paths=[], **kwargs):
    ctx = create_test_context()
    with TempDirectoryContext() as tdc:
        test_data = os.path.join(tdc.temp_directory, "test-data")
        os.makedirs(test_data)
        kwargs["test_data"] = test_data
        with galaxy_config(ctx, tool_paths, **kwargs) as gc:
            yield gc
