"""Tests for the ``planemo.galaxy.serve`` module.

This tests this as a library functionality - additional integration
style tests are available in ``test_cmd_serve.py``.
"""
import os

from planemo import network_util
from planemo import shed
from planemo.galaxy import galaxy_serve, shed_serve
from planemo.runnable import for_path

from .test_utils import (
    CliTestCase,
    PROJECT_TEMPLATES_DIR,
    skip_if_environ,
    TEST_DATA_DIR,
    TEST_REPOS_DIR,
)


class GalaxyServeTestCase(CliTestCase):
    """Tests for planemo.galaxy.serve."""

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    def test_serve_daemon(self):
        """Test serving a galaxy tool via a daemon Galaxy process."""
        port = network_util.get_free_port()
        cat_path = os.path.join(TEST_REPOS_DIR, "single_tool", "cat.xml")
        config = galaxy_serve(
            self.test_context,
            [for_path(cat_path)],
            install_galaxy=True,
            port=port,
            daemon=True,
        )

        assert network_util.wait_net_service(
            "localhost",
            config.port,
            timeout=.1,
        )
        config_dict = config.gi.config.get_config()
        assert "allow_user_dataset_purge" in config_dict
        config.kill()
        assert not network_util.wait_net_service(
            "localhost",
            config.port,
            timeout=.1,
        )

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    def test_serve_workflow(self):
        """Test serving a galaxy workflow via a daemon Galaxy process."""
        port = network_util.get_free_port()
        random_lines = os.path.join(PROJECT_TEMPLATES_DIR, "demo", "randomlines.xml")
        cat = os.path.join(PROJECT_TEMPLATES_DIR, "demo", "cat.xml")
        worklfow = os.path.join(TEST_DATA_DIR, "wf1.gxwf.yml")
        extra_tools = [random_lines, cat]
        config = galaxy_serve(
            self.test_context,
            [for_path(worklfow)],
            install_galaxy=True,
            port=port,
            daemon=True,
            extra_tools=extra_tools,
        )
        user_gi = config.user_gi
        assert user_gi.tools.get_tools(tool_id="random_lines1")
        assert len(user_gi.workflows.get_workflows()) == 1
        config.kill()

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    def test_shed_serve_daemon(self):
        """Test serving FASTQC from the tool shed via a daemon Galaxy process."""
        port = network_util.get_free_port()
        fastqc_path = os.path.join(TEST_REPOS_DIR, "fastqc")
        ctx = self.test_context
        install_args_list = shed.install_arg_lists(
            ctx, [fastqc_path],
            shed_target="toolshed",
        )
        with shed_serve(
            ctx, install_args_list,
            port=port,
            skip_dependencies=True,
            install_galaxy=True,
        ) as config:
            assert network_util.wait_net_service(
                "localhost",
                config.port,
                timeout=.1,
            )
