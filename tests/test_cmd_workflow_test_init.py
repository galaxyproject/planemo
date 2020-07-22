"""Tests for the ``workflow_test_init`` command."""
import os

import yaml

from .test_utils import (
    CliTestCase,
)


class CmdWorkflowTestInitTestCase(CliTestCase):

    def test_plain_init(self):
        with self._isolate_with_test_data("wf_repos/from_format2/0_basic_native") as f:
            init_cmd = ["workflow_test_init", "0_basic_native.yml"]
            self._check_exit_code(init_cmd)
            test_path = os.path.join(f, "0_basic_native_tests.yml")
            assert os.path.exists(test_path)
            with open(test_path, "r") as stream:
                test_config = yaml.safe_load(stream)
            assert isinstance(test_config, list)
            job_path = os.path.join(f, "0_basic_native_job1.yml")
            assert os.path.exists(job_path)
            with open(job_path, "r") as stream:
                job = yaml.safe_load(stream)
            assert isinstance(job, dict)
