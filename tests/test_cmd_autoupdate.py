"""Tests for the ``autoupdate`` command."""
import tempfile

from .test_utils import (
    CliTestCase
)


def create_tmp_test_tool_file():
    xml_str = b"""<tool id="autoupdate_test" version="@TOOL_VERSION@+galaxy@VERSION_SUFFIX@">
    <macros>
        <token name="@TOOL_VERSION@">0.6.0</token>
        <token name="@VERSION_SUFFIX@">1</token>
    </macros>
    <requirements>
        <requirement type="package" version="@TOOL_VERSION@">xopen</requirement>
        <requirement type="package" version="3.7.0">zeroc-ice</requirement>
    </requirements>
</tool>
    """
    t = tempfile.NamedTemporaryFile(suffix='.xml', delete=False)
    t.write(xml_str)
    return t.name


class CmdAutoupdateTestCase(CliTestCase):
    """Container class defining test cases for the ``autoupdate`` command."""

    def test_autoupdate_dry_run(self):
        """Test autoupdate command with dry run flag."""
        with self._isolate():
            xmlfile = create_tmp_test_tool_file()
            autoupdate_command = [
                "autoupdate",
                xmlfile,
                "--conda_channels", "bioconda",
                "--dry-run"
            ]
            self._check_exit_code(autoupdate_command, exit_code=0)

    def test_autoupdate(self):
        """Test autoupdate command."""
        with self._isolate():
            xmlfile = create_tmp_test_tool_file()
            autoupdate_command = [
                "autoupdate",
                xmlfile,
                "--conda_channels", "bioconda"
            ]
            self._check_exit_code(autoupdate_command, exit_code=0)
            with open(xmlfile) as stream:
                updated_tool = stream.readlines()
                assert updated_tool[2].strip() == '<token name="@TOOL_VERSION@">0.7.3</token>'
                assert updated_tool[3].strip() == '<token name="@VERSION_SUFFIX@">0</token>'
                assert updated_tool[7].strip() == '<requirement type="package" version="3.7.1">zeroc-ice</requirement>'
