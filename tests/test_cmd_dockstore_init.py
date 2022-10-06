"""Tests for the ``dockstore_init`` command."""
import os
from typing import Optional

import yaml

from .test_utils import CliTestCase


class CmdDockstoreInitTestCase(CliTestCase):
    def test_plain_init(self):
        self.run_dockstore_init()

    def test_init_publish_true(self):
        self.run_dockstore_init(True)

    def test_init_publish_false(self):
        self.run_dockstore_init(False)

    def run_dockstore_init(self, publish: Optional[bool] = None):
        with self._isolate_with_test_data("wf_repos/from_format2/0_basic_native") as f:
            init_cmd = ["dockstore_init"]
            expect_published = True
            if publish:
                init_cmd.append("--publish")
            elif publish is False:
                expect_published = False
                init_cmd.append("--no_publish")
            self._check_exit_code(init_cmd)
            assert os.path.exists(".dockstore.yml")
            with open(".dockstore.yml") as fh:
                dockstore_config = yaml.safe_load(fh)
            assert str(dockstore_config["version"]) == "1.2"
            assert "workflows" in dockstore_config
            assert len(dockstore_config["workflows"]) == 1
            assert dockstore_config["workflows"][0]["publish"] == expect_published
            workflow_lint_cmd = ["workflow_lint", "--fail_level", "error", f]
            self._check_exit_code(workflow_lint_cmd)
