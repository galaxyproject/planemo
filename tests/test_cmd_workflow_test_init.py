"""Tests for the ``workflow_test_init`` command."""
import os

import yaml

from .test_utils import CliTestCase


class CmdWorkflowTestInitTestCase(CliTestCase):
    def test_plain_init(self):
        with self._isolate_with_test_data("wf_repos/from_format2/0_basic_native") as f:
            init_cmd = ["workflow_test_init", "0_basic_native.yml"]
            self._check_exit_code(init_cmd)
            test_path = os.path.join(f, "0_basic_native_tests.yml")
            assert os.path.exists(test_path)
            with open(test_path) as stream:
                test_config = yaml.safe_load(stream)
            job = test_config[0]["job"]
            assert isinstance(job, dict)

    def test_init_split_test(self):
        with self._isolate_with_test_data("wf_repos/from_format2/0_basic_native") as f:
            init_cmd = ["workflow_test_init", "--split_test", "0_basic_native.yml"]
            self._check_exit_code(init_cmd)
            test_path = os.path.join(f, "0_basic_native_tests.yml")
            assert os.path.exists(test_path)
            with open(test_path) as stream:
                test_config = yaml.safe_load(stream)
            assert isinstance(test_config, list)
            job_path = os.path.join(f, "0_basic_native_job1.yml")
            assert os.path.exists(job_path)
            with open(job_path) as stream:
                job = yaml.safe_load(stream)
            assert isinstance(job, dict)
            assert "the_input" in job
            assert job.get("the_input").get("path") == "todo_test_data_path.ext"

    def test_cannot_overwrite(self):
        with self._isolate_with_test_data("wf_repos/from_format2/0_basic_native") as f:
            init_cmd = ["workflow_test_init", "0_basic_native.yml"]
            test_path = os.path.join(f, "0_basic_native_tests.yml")
            with open(test_path, "w") as f:
                f.write("already exists")
            self._check_exit_code(init_cmd, exit_code=1)

    def test_force(self):
        with self._isolate_with_test_data("wf_repos/from_format2/0_basic_native") as f:
            init_cmd = ["workflow_test_init", "--force", "0_basic_native.yml"]
            test_path = os.path.join(f, "0_basic_native_tests.yml")
            with open(test_path, "w") as f:
                f.write("already exists")
            self._check_exit_code(init_cmd)
