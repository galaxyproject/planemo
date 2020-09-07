"""Tests for the ``autoupdate`` command."""
import os
from shutil import copyfile

from .test_utils import (
    CliTestCase,
    TEST_DATA_DIR
)


class CmdAutoupdateTestCase(CliTestCase):
    """Container class defining test cases for the ``autoupdate`` command."""

    # @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    def test_autoupdate_dry_run(self):
        """Test autoupdate command with dry run flag."""
        with self._isolate():
            training_init_command = [
                "autoupdate",
                "{}/autoupdate_test.xml".format(TEST_DATA_DIR),
                "--conda_channels", "bioconda",
                "--dry-run"
            ]
            self._check_exit_code(training_init_command, exit_code=0)

    def test_autoupdate(self):
        """Test autoupdate command."""
        copyfile("{}/autoupdate_test.xml".format(TEST_DATA_DIR),  "{}/autoupdate_test_tmp.xml".format(TEST_DATA_DIR))
        with self._isolate():
            training_init_command = [
                "autoupdate",
                "{}/autoupdate_test_tmp.xml".format(TEST_DATA_DIR),
                "--conda_channels", "bioconda"
            ]

            self._check_exit_code(training_init_command, exit_code=0)

        with open("{}/autoupdate_test_tmp.xml".format(TEST_DATA_DIR)) as f:
            updated_tool = f.readlines()
            assert updated_tool[2].strip() == '<token name="@TOOL_VERSION@">0.7.3</token>'
            assert updated_tool[3].strip() == '<token name="@GALAXY_VERSION@">0</token>'
            assert updated_tool[7].strip() == '<requirement type="package" version="3.7.1">zeroc-ice</requirement>'
        os.remove("{}/autoupdate_test_tmp.xml".format(TEST_DATA_DIR))
