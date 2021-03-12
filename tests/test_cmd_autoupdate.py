"""Tests for the ``autoupdate`` command."""
import tempfile

from .test_utils import (
    CliTestCase
)


def create_tmp_test_tool_file(tool_version):
    """
    Note that to ensure this test is stable, we use packages that have
    been migrated from bioconda to conda-forge and test with --conda_channels bioconda
    so that the versions are guaranteed not to increase in future.
    """
    xml_str = """<tool id="autoupdate_test" version="@TOOL_VERSION@+galaxy@VERSION_SUFFIX@">
    <macros>
        <token name="@TOOL_VERSION@">{}</token>
        <token name="@VERSION_SUFFIX@">1</token>
    </macros>
    <requirements>
        <requirement type="package" version="@TOOL_VERSION@">xopen</requirement>
        <requirement type="package" version="2015">smina</requirement>
    </requirements>
</tool>
    """.format(tool_version)
    t = tempfile.NamedTemporaryFile(suffix='.xml', delete=False, mode='w')
    t.write(xml_str)
    return t.name


class CmdAutoupdateTestCase(CliTestCase):
    """Container class defining test cases for the ``autoupdate`` command."""

    def setUp(self):
        super(CmdAutoupdateTestCase, self).setUp()
        self._runner.invoke(self._cli.planemo, ['conda_init'])

    def test_autoupdate_dry_run(self):
        """Test autoupdate command with dry run flag."""
        with self._isolate():
            xmlfile = create_tmp_test_tool_file('0.6.0')
            autoupdate_command = [
                "autoupdate",
                xmlfile,
                "--conda_channels", "bioconda",
                "--dry-run"
            ]
            result = self._runner.invoke(self._cli.planemo, autoupdate_command)
            assert "Update required to {}!".format(xmlfile) in result.output
            assert "Tool main requirement has version 0.6.0, newest conda version is 0.7.3" in result.output

    def test_autoupdate(self):
        """Test autoupdate command."""
        with self._isolate():
            xmlfile = create_tmp_test_tool_file('0.6.0')
            autoupdate_command = [
                "autoupdate",
                xmlfile,
                "--conda_channels", "bioconda"
            ]
            result = self._runner.invoke(self._cli.planemo, autoupdate_command)
            assert 'Tool {} updated.'.format(xmlfile) in result.output
            with open(xmlfile) as f:
                xmlfile_contents = f.read()
            assert "2017.11.9" in xmlfile_contents

    def test_autoupdate_no_update_needed(self):
        """Test autoupdate command when no update is needed."""
        with self._isolate():
            xmlfile = create_tmp_test_tool_file('0.7.3')
            autoupdate_command = [
                "autoupdate",
                xmlfile,
                "--conda_channels", "bioconda"
            ]
            result = self._runner.invoke(self._cli.planemo, autoupdate_command)
            assert 'No updates required or made to {}.'.format(xmlfile) in result.output
