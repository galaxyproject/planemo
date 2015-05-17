import os

from .test_utils import (
    CliTestCase,
    skip_if_environ,
    TEST_REPOS_DIR,
)
from . import network_util
from planemo import galaxy_serve


class GalaxyServeTestCase(CliTestCase):

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    def test_serve_daemon(self):
        port = network_util.get_free_port()
        cat_path = os.path.join(TEST_REPOS_DIR, "single_tool", "cat.xml")
        config = galaxy_serve.serve(
            self.test_context,
            cat_path,
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
