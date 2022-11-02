"""Tests for the ``autoupdate`` command."""
import json
import os
import tempfile
from contextlib import contextmanager

import yaml

from .test_utils import CliTestCase


@contextmanager
def create_tmp_test_tool_file(tool_version):
    """
    Note that to ensure this test is stable, we use packages that have
    been migrated from bioconda to conda-forge and test with --conda_channels bioconda
    so that the versions are guaranteed not to increase in future.
    """
    xml_str = f"""<tool id="autoupdate_test" version="@TOOL_VERSION@+galaxy@VERSION_SUFFIX@">
    <macros>
        <token name="@TOOL_VERSION@">{tool_version}</token>
        <token name="@VERSION_SUFFIX@">1</token>
    </macros>
    <requirements>
        <requirement type="package" version="@TOOL_VERSION@">xopen</requirement>
        <requirement type="package" version="2015">smina</requirement>
    </requirements>
</tool>
    """
    with tempfile.TemporaryDirectory() as tempdir, tempfile.NamedTemporaryFile(
        suffix=".xml", mode="w", dir=tempdir
    ) as t:
        t.write(xml_str)
        t.flush()
        yield t.name


class CmdAutoupdateTestCase(CliTestCase):
    """Container class defining test cases for the ``autoupdate`` command."""

    def setUp(self):
        super().setUp()

    def test_autoupdate_dry_run(self):
        """Test autoupdate command with dry run flag."""
        with self._isolate(), create_tmp_test_tool_file("0.6.0") as xmlfile:
            autoupdate_command = ["autoupdate", xmlfile, "--conda_channels", "bioconda", "--dry-run"]
            result = self._runner.invoke(self._cli.planemo, autoupdate_command)
            assert f"Update required to {xmlfile}!" in result.output
            assert "Tool main requirement has version 0.6.0, newest conda version is 0.7.3" in result.output

    def test_autoupdate(self):
        """Test autoupdate command."""
        with self._isolate(), create_tmp_test_tool_file("0.6.0") as xmlfile:
            autoupdate_command = ["autoupdate", xmlfile, "--conda_channels", "bioconda"]
            result = self._runner.invoke(self._cli.planemo, autoupdate_command)
            assert f'Updating {xmlfile.split("/")[-1]} from version 0.6.0 to 0.7.3' in result.output
            assert f"Tool {xmlfile} successfully updated." in result.output
            with open(xmlfile) as f:
                xmlfile_contents = f.read()
            assert "2017.11.9" in xmlfile_contents

    def test_autoupdate_directory(self):
        """Test autoupdate command."""
        with self._isolate(), create_tmp_test_tool_file("0.6.0") as xmlfile:
            xml_directory = os.path.dirname(xmlfile)
            autoupdate_command = ["autoupdate", xml_directory, "--conda_channels", "bioconda"]
            result = self._runner.invoke(self._cli.planemo, autoupdate_command)
            assert f'Updating {xmlfile.split("/")[-1]} from version 0.6.0 to 0.7.3' in result.output
            assert f"Tool {xmlfile} successfully updated." in result.output
            with open(xmlfile) as f:
                xmlfile_contents = f.read()
            assert "2017.11.9" in xmlfile_contents

    def test_autoupdate_no_update_needed(self):
        """Test autoupdate command when no update is needed."""
        with self._isolate(), create_tmp_test_tool_file("0.7.3") as xmlfile:
            autoupdate_command = ["autoupdate", xmlfile, "--conda_channels", "bioconda"]
            result = self._runner.invoke(self._cli.planemo, autoupdate_command)
            assert f"No updates required or made to {xmlfile}." in result.output

    def test_autoupdate_workflow(self):
        """Test autoupdate command for a workflow is needed."""
        with self._isolate_with_test_data("wf_repos/autoupdate_tests") as f:
            wf_file = os.path.join(f, "diff-refactor-test.ga")
            autoupdate_command = ["autoupdate", wf_file, "--galaxy_branch", "dev"]  # need >= 21.05
            result = self._runner.invoke(self._cli.planemo, autoupdate_command)

            assert f"Auto-updating workflow {wf_file}" in result.output
            with open(wf_file) as g:
                wf = json.load(g)
            # check tool within parent wf has updated
            assert wf["steps"]["1"]["tool_version"] == "3.7+galaxy0"
            # check tool within subworkflow has updated
            assert wf["steps"]["2"]["subworkflow"]["steps"]["1"]["tool_version"] == "3.7+galaxy0"
            assert (
                wf["steps"]["2"]["subworkflow"]["steps"]["1"]["tool_id"]
                == "toolshed.g2.bx.psu.edu/repos/bgruening/diff/diff/3.7+galaxy0"
            )
            assert wf["version"] == 2
            assert wf["release"] == "0.1.1"

            result = self._runner.invoke(self._cli.planemo, autoupdate_command)  # rerun on already updated WF
            assert "No newer tool versions were found, so the workflow was not updated." in result.output

            wf_file = os.path.join(f, "diff-refactor-test.gxwf.yml")
            autoupdate_command[1] = wf_file
            result = self._runner.invoke(self._cli.planemo, autoupdate_command)
            assert f"Auto-updating workflow {wf_file}" in result.output

            with open(wf_file) as f:
                wf = yaml.safe_load(f)
            assert wf["steps"][0]["tool_version"] == "3.7+galaxy0"
            assert wf["steps"][0]["tool_id"] == "toolshed.g2.bx.psu.edu/repos/bgruening/diff/diff/3.7+galaxy0"
            assert wf["steps"][1]["run"]["steps"][0]["tool_version"] == "3.7+galaxy0"
            assert wf["release"] == "0.1.1"

    def test_autoupdate_workflow_unexisting_version(self):
        """Test autoupdate command for a workflow where the version of the tool is not in the toolshed."""
        with self._isolate_with_test_data("wf_repos/autoupdate_tests") as f:
            wf_file = os.path.join(f, "workflow_with_unexisting_version_of_tool.ga")
            autoupdate_command = ["autoupdate", wf_file]
            self._runner.invoke(self._cli.planemo, autoupdate_command)
            # We just want to be sure planemo autoupdate do not raise an error
            # Currently it would write to the output that no update are available
            # In future versions it could be great that it gives the last valid version.

    def test_autoupdate_workflow_unexisting_tool(self):
        """Test autoupdate command for a workflow where the tool is not in the toolshed."""
        with self._isolate_with_test_data("wf_repos/autoupdate_tests") as f:
            wf_file = os.path.join(f, "workflow_with_unexisting_tool.ga")
            autoupdate_command = ["autoupdate", wf_file]
            result = self._runner.invoke(self._cli.planemo, autoupdate_command)
            assert "No newer tool versions were found, so the workflow was not updated." in result.output
