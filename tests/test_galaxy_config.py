"""Unit tests for ``planemo.galaxy.config``."""

import contextlib
import json
import os

import yaml

from planemo.galaxy.config import (
    galaxy_config,
    get_refgenie_config,
    write_galaxy_config,
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


def test_gxits_enabled_by_default():
    """Test that Galaxy Interactive Tools are enabled by default in galaxy.yml."""
    with _test_write_galaxy_config() as (config_data, properties, env):
        gravity = config_data["gravity"]
        gx_it_proxy = gravity["gx_it_proxy"]
        assert gx_it_proxy["enable"] is True
        assert "port" in gx_it_proxy
        assert isinstance(gx_it_proxy["port"], int)


def test_gxits_sets_galaxy_properties():
    """Test that enabling GxITs sets the required Galaxy properties."""
    with _test_write_galaxy_config() as (config_data, properties, env):
        assert properties["interactivetools_enable"] is True
        assert properties["interactivetools_upstream_proxy"] is False
        assert "galaxy_infrastructure_url" in properties
        assert "interactivetools_proxy_host" in properties
        # The proxy host should include the gx_it_proxy port
        gx_it_port = config_data["gravity"]["gx_it_proxy"]["port"]
        assert properties["interactivetools_proxy_host"] == f"localhost:{gx_it_port}"


def test_gxits_disabled_with_flag():
    """Test that --disable_gxits flag properly disables interactive tools."""
    with _test_write_galaxy_config(disable_gxits=True) as (config_data, properties, env):
        gravity = config_data["gravity"]
        gx_it_proxy = gravity["gx_it_proxy"]
        assert gx_it_proxy["enable"] is False
        assert "port" not in gx_it_proxy
        # Galaxy properties for interactive tools should not be set
        assert "interactivetools_enable" not in properties
        assert "interactivetools_upstream_proxy" not in properties
        assert "interactivetools_proxy_host" not in properties


def test_gxits_infrastructure_url_uses_host_and_port():
    """Test that the galaxy_infrastructure_url uses the configured host and port."""
    with _test_write_galaxy_config(host="0.0.0.0", port=9999) as (config_data, properties, env):
        assert properties["galaxy_infrastructure_url"] == "http://0.0.0.0:9999"
        assert properties["interactivetools_proxy_host"].startswith("0.0.0.0:")


@contextlib.contextmanager
def _test_galaxy_config(tool_paths=[], **kwargs):
    ctx = create_test_context()
    with TempDirectoryContext() as tdc:
        test_data = os.path.join(tdc.temp_directory, "test-data")
        os.makedirs(test_data)
        kwargs["test_data"] = test_data
        with galaxy_config(ctx, tool_paths, **kwargs) as gc:
            yield gc


@contextlib.contextmanager
def _test_write_galaxy_config(**kwds):
    """Helper to test write_galaxy_config directly with a fake galaxy root."""
    with TempDirectoryContext() as tdc:
        galaxy_root = tdc.temp_directory
        # Create fake galaxy version file (>= 22.01 to use YAML config)
        galaxy_lib_path = os.path.join(galaxy_root, "lib", "galaxy")
        os.makedirs(galaxy_lib_path)
        version_path = os.path.join(galaxy_lib_path, "version.py")
        with open(version_path, "w") as version_fh:
            version_fh.write('VERSION_MAJOR = "24.1"')

        config_dir = os.path.join(galaxy_root, "config")
        os.makedirs(config_dir)

        def config_join(*args):
            return os.path.join(config_dir, *args)

        host = kwds.pop("host", "localhost")
        port = kwds.pop("port", 9090)
        properties = {}
        env = {}
        template_args = {"port": port}
        config_kwds = {"host": host}
        config_kwds.update(kwds)

        write_galaxy_config(
            galaxy_root=galaxy_root,
            properties=properties,
            env=env,
            kwds=config_kwds,
            template_args=template_args,
            config_join=config_join,
        )

        config_file = env["GALAXY_CONFIG_FILE"]
        with open(config_file) as f:
            config_data = json.load(f)

        yield config_data, properties, env
