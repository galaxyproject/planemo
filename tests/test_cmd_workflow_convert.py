"""Tests for the ``workflow_test_init`` command."""
import json
import os

import yaml

from .test_utils import (
    CliTestCase,
    TEST_DATA_DIR,
)


class CmdWorkflowConvertTestCase(CliTestCase):
    def test_gxwf_to_ga(self):
        with self._isolate() as f:
            gx2_wf_path = os.path.join(TEST_DATA_DIR, "wf1.gxwf.yml")
            init_cmd = ["workflow_convert", gx2_wf_path, "-o", "wf1-test.ga"]
            self._check_exit_code(init_cmd)
            ga_wf_path = os.path.join(f, "wf1-test.ga")

            assert os.path.exists(ga_wf_path)
            with open(ga_wf_path) as stream:
                wf = json.load(stream)
            assert isinstance(wf, dict)
            assert wf["steps"]["1"]["tool_id"] == "cat"

    def test_ga_to_gxwf(self):
        with self._isolate() as f:
            ga_wf_path = os.path.join(TEST_DATA_DIR, "wf1.ga")
            init_cmd = ["workflow_convert", ga_wf_path, "-o", "wf1-test.yml"]
            self._check_exit_code(init_cmd)
            gx2_wf_path = os.path.join(f, "wf1-test.yml")

            assert os.path.exists(gx2_wf_path)
            with open(gx2_wf_path) as stream:
                wf = yaml.safe_load(stream)
            assert isinstance(wf, dict)
            assert wf["steps"]["first_cat"]["tool_id"] == "cat"
