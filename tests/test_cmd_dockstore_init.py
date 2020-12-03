"""Tests for the ``dockstore_init`` command."""
import os

import yaml

from .test_utils import (
    CliTestCase,
)


class CmdDockstoreInitTestCase(CliTestCase):

    def test_plain_init(self):
        with self._isolate_with_test_data("wf_repos/from_format2/0_basic_native") as f:
            init_cmd = ["dockstore_init"]
            self._check_exit_code(init_cmd)
            assert os.path.exists(".dockstore.yml")
            with open(".dockstore.yml", "r") as f:
                dockstore_config = yaml.safe_load(f)
            assert str(dockstore_config["version"]) == "1.2"
            assert "workflows" in dockstore_config
            assert len(dockstore_config["workflows"]) == 1
