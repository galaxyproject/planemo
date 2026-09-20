"""Tests for the ``profile_job_config_init`` command."""

import json
import os
from tempfile import mkdtemp

import yaml

from .test_utils import CliTestCase

TEST_PROFILE_NAME = "profilejobtest"


class CmdJobConfigInitTestCase(CliTestCase):
    def test_job_config_init_simple(self):
        with self._isolate():
            self._check_exit_code(["profile_create", TEST_PROFILE_NAME])
            init_cmd = ["profile_job_config_init", TEST_PROFILE_NAME]
            self._check_exit_code(init_cmd)
            self._check_exit_code(["profile_delete", TEST_PROFILE_NAME])

    def test_job_config_recorded_in_profile(self):
        with self._isolate():
            workspace = mkdtemp()
            self._check_exit_code(["--directory", workspace, "profile_create", TEST_PROFILE_NAME])
            self._check_exit_code(
                ["--directory", workspace, "profile_job_config_init", TEST_PROFILE_NAME, "--runner", "slurm"]
            )
            profile_directory = os.path.join(workspace, "profiles", TEST_PROFILE_NAME)
            with open(os.path.join(profile_directory, "planemo_profile_options.json")) as fh:
                profile_options = json.load(fh)

            job_config_file = profile_options["job_config_file"]
            assert job_config_file == os.path.join(profile_directory, "job_conf.yml")
            with open(job_config_file) as fh:
                job_config = yaml.safe_load(fh)
            assert job_config["execution"]["default"] == "slurm"
