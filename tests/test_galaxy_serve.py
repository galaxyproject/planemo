"""Tests for the ``planemo.galaxy.serve`` module.

This tests this as a library functionality - additional integration
style tests are available in ``test_cmd_serve.py``.
"""
import os

from planemo import network_util
from planemo.galaxy import galaxy_serve
from planemo.runnable import for_path
from .test_utils import (
    CliTestCase,
    mark,
    PROJECT_TEMPLATES_DIR,
    skip_if_environ,
    target_galaxy_branch,
    TEST_DATA_DIR,
    TEST_REPOS_DIR,
)


class GalaxyServeTestCase(CliTestCase):
    """Tests for planemo.galaxy.serve."""

    @skip_if_environ("PLANEMO_SKIP_REDUNDANT_TESTS")  # redundant with test_cmd_serve -> test_serve_daemon
    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    @mark.tests_galaxy_branch
    def test_serve_daemon(self):
        """Test serving a galaxy tool via a daemon Galaxy process."""
        port = network_util.get_free_port()
        cat_path = os.path.join(TEST_REPOS_DIR, "single_tool", "cat.xml")
        with galaxy_serve(
            self.test_context,
            [for_path(cat_path)],
            install_galaxy=True,
            galaxy_branch=target_galaxy_branch(),
            port=port,
            daemon=True,
            no_dependency_resolution=True,
        ) as config:
            _assert_service_up(config)
            config.kill()
        _assert_service_down(config)

    @skip_if_environ("PLANEMO_SKIP_REDUNDANT_TESTS")  # redundant with test_cmd_serve -> test_serve_workflow
    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    @mark.tests_galaxy_branch
    def test_serve_workflow(self):
        """Test serving a galaxy workflow via a daemon Galaxy process."""
        port = network_util.get_free_port()
        random_lines = os.path.join(PROJECT_TEMPLATES_DIR, "demo", "randomlines.xml")
        cat = os.path.join(PROJECT_TEMPLATES_DIR, "demo", "cat.xml")
        workflow = os.path.join(TEST_DATA_DIR, "wf1.gxwf.yml")
        extra_tools = [random_lines, cat]
        with galaxy_serve(
            self.test_context,
            [for_path(workflow)],
            install_galaxy=True,
            galaxy_branch=target_galaxy_branch(),
            port=port,
            daemon=True,
            extra_tools=extra_tools,
            no_dependency_resolution=True,
        ) as config:
            _assert_service_up(config)
            user_gi = config.user_gi
            assert user_gi.tools.show_tool("random_lines1")
            assert len(user_gi.workflows.get_workflows()) == 1
            config.kill()
            _assert_service_down(config)


def _assert_service_up(config):
    assert network_util.wait_net_service(
        "localhost",
        config.port,
        timeout=0.1,
    )
    galaxy_config_api = config.gi.config
    config_dict = galaxy_config_api.get_config()
    assert "allow_user_dataset_purge" in config_dict


def _assert_service_down(config):
    assert not network_util.wait_net_service(
        "localhost",
        config.port,
        timeout=0.1,
    )
