"""Unit tests for ``planemo.galaxy.config``."""
import contextlib
import os

from .test_utils import TempDirectoryContext, test_context

from planemo.galaxy.config import galaxy_config


def test_defaults():
    """Test by default Galaxy files are stored in temp ``config_directory``."""
    with _test_galaxy_config() as config:
        config_directory = config.config_directory
        _assert_property_is(config, "file_path", os.path.join(config_directory, "files"))
        conn = "sqlite:///%s/galaxy.sqlite?isolation_level=IMMEDIATE" % config_directory
        _assert_property_is(config, "database_connection", conn)


def test_database_connection_override_path():
    """Test by default Galaxy files are stored in temp ``config_directory``."""
    conn = "postgresql://username:password@localhost/mydatabase"
    with _test_galaxy_config(database_connection=conn) as config:
        _assert_property_is(config, "database_connection", conn)


def test_override_files_path():
    """Test Galaxy file path overrideable with --file_path."""
    with TempDirectoryContext() as tdc:
        with _test_galaxy_config(file_path=tdc.temp_directory) as config:
            _assert_property_is(config, "file_path", tdc.temp_directory)


def _assert_property_is(config, prop, value):
    env_var = "GALAXY_CONFIG_OVERRIDE_%s" % prop.upper()
    assert config.env[env_var] == value


@contextlib.contextmanager
def _test_galaxy_config(tool_paths=[], **kwargs):
    ctx = test_context()
    with TempDirectoryContext() as tdc:
        test_data = os.path.join(tdc.temp_directory, "test-data")
        os.makedirs(test_data)
        kwargs["test_data"] = test_data
        with galaxy_config(ctx, tool_paths, **kwargs) as gc:
            yield gc
