"""Tests for the ``job_config_init`` command."""

from .test_utils import CliTestCase

TEST_PROFILE_NAME = "profilejobtest"


class CmdJobConfigInitTestCase(CliTestCase):
    def test_job_config_init_simple(self):
        with self._isolate():
            self._check_exit_code(["profile_create", TEST_PROFILE_NAME])
            init_cmd = ["profile_job_config_init", TEST_PROFILE_NAME]
            self._check_exit_code(init_cmd)
            self._check_exit_code(["profile_delete", TEST_PROFILE_NAME])
