"""Tests for the ``autoupdate`` command."""

import json
import os
import shutil
import tempfile
from contextlib import contextmanager

import yaml

from .test_utils import (
    CliTestCase,
    skip_if_environ,
)


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
        yield os.path.realpath(t.name)


class CmdAutoupdateTestCase(CliTestCase):
    """Container class defining test cases for the ``autoupdate`` command."""

    def setUp(self):
        super().setUp()

    @skip_if_environ("PLANEMO_SKIP_SLOW_TESTS")
    def test_autoupdate_dry_run(self):
        """Test autoupdate command with dry run flag."""
        with self._isolate(), create_tmp_test_tool_file("0.6.0") as xmlfile:
            autoupdate_command = ["autoupdate", xmlfile, "--conda_channels", "bioconda", "--dry-run"]
            result = self._runner.invoke(self._cli.planemo, autoupdate_command)
            assert f"Update required to {xmlfile}!" in result.output
            assert "Tool main requirement has version 0.6.0, newest conda version is 0.7.3" in result.output

    @skip_if_environ("PLANEMO_SKIP_SLOW_TESTS")
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

    @skip_if_environ("PLANEMO_SKIP_SLOW_TESTS")
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

    @skip_if_environ("PLANEMO_SKIP_SLOW_TESTS")
    def test_autoupdate_no_update_needed(self):
        """Test autoupdate command when no update is needed."""
        with self._isolate(), create_tmp_test_tool_file("0.7.3") as xmlfile:
            autoupdate_command = ["autoupdate", xmlfile, "--conda_channels", "bioconda"]
            result = self._runner.invoke(self._cli.planemo, autoupdate_command)
            assert f"No updates required or made to {xmlfile}." in result.output

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    def test_autoupdate_multiple_workflows(self):
        """Test autoupdate command for a workflow is needed."""
        with self._isolate_with_test_data("wf_repos/autoupdate_tests") as f, tempfile.TemporaryDirectory(
            dir=f, prefix="autoupdate_test"
        ) as isolated_dir:
            source_file = os.path.join(f, "diff-refactor-test.ga")
            # We update identical workflows in the same autoupdate call,
            # both workflows must be updated.
            targets = [os.path.join(isolated_dir, wf) for wf in ("wf1.ga", "wf2.ga")]
            for target in targets:
                shutil.copy(source_file, target)
            autoupdate_command = ["autoupdate", *targets]
            result = self._runner.invoke(self._cli.planemo, autoupdate_command)
            assert "Auto-updating workflow" in result.output
            for wf_file in targets:
                with open(wf_file) as g:
                    wf = json.load(g)
                # check tool within parent wf has updated
                tool_step = next(iter(step for step in wf["steps"].values() if step["type"] == "tool"))
                assert tool_step["tool_version"] != "3.6+galaxy1"
                subworkflow_step = next(iter(step for step in wf["steps"].values() if step["type"] == "subworkflow"))
                # check tool within subworkflow has updated
                assert subworkflow_step["subworkflow"]["steps"]["1"]["tool_version"] != "3.6+galaxy1"
                assert (
                    subworkflow_step["subworkflow"]["steps"]["1"]["tool_id"]
                    != "toolshed.g2.bx.psu.edu/repos/bgruening/diff/diff/3.6+galaxy1"
                )
                assert wf["release"] == "0.1.1"

                result = self._runner.invoke(self._cli.planemo, autoupdate_command)  # rerun on already updated WF
                assert "No newer tool versions were found, so the workflow was not updated." in result.output

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    def test_autoupdate_gxformat2_workflow(self):
        with self._isolate_with_test_data("wf_repos/autoupdate_tests") as f:
            wf_file = os.path.join(f, "diff-refactor-test.gxwf.yml")
            autoupdate_command = ["autoupdate", wf_file]
            result = self._runner.invoke(self._cli.planemo, autoupdate_command)
            assert f"Auto-updating workflow {wf_file}" in result.output

            with open(wf_file) as f:
                wf = yaml.safe_load(f)

            tool_step = next(iter(step for step in wf["steps"] if "tool_id" in step))
            assert tool_step["tool_version"] != "3.6+galaxy1"
            assert tool_step["tool_id"] != "toolshed.g2.bx.psu.edu/repos/bgruening/diff/diff/3.6+galaxy1"
            workflow_step = next(iter(step for step in wf["steps"] if "run" in step))
            assert workflow_step["run"]["steps"][0]["tool_version"] != "3.6+galaxy1"
            assert wf["release"] == "0.1.1"

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    def test_autoupdate_workflow_from_multiple_tool_sheds(self):
        with self._isolate_with_test_data("wf_repos/autoupdate_tests") as f:
            wf_file = os.path.join(f, "wf_autoupdate_test_multiple_repos.ga")
            autoupdate_command = ["autoupdate", wf_file]
            result = self._runner.invoke(self._cli.planemo, autoupdate_command)
            assert f"Auto-updating workflow {wf_file}" in result.output
            with open(wf_file) as g:
                wf = json.load(g)
            # Assert toolshed tool is updated
            assert wf["steps"]["1"]["tool_version"] != "9.3+galaxy0"
            # Assert testtoolshed tool is updated
            assert wf["steps"]["2"]["tool_version"] != "0.69"

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    def test_autoupdate_workflow_unexisting_version(self):
        """Test autoupdate command for a workflow where the version of the tool is not in the toolshed."""
        with self._isolate_with_test_data("wf_repos/autoupdate_tests") as f:
            wf_file = os.path.join(f, "workflow_with_unexisting_version_of_tool.ga")
            autoupdate_command = ["autoupdate", wf_file]
            self._runner.invoke(self._cli.planemo, autoupdate_command)
            # We just want to be sure planemo autoupdate do not raise an error
            # Currently it would write to the output that no update are available
            # In future versions it could be great that it gives the last valid version.

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    def test_autoupdate_workflow_unexisting_tool(self):
        """Test autoupdate command for a workflow where the tool is not in the toolshed."""
        with self._isolate_with_test_data("wf_repos/autoupdate_tests") as f:
            wf_file = os.path.join(f, "workflow_with_unexisting_tool.ga")
            autoupdate_command = ["autoupdate", wf_file]
            result = self._runner.invoke(self._cli.planemo, autoupdate_command)
            assert "No newer tool versions were found, so the workflow was not updated." in result.output
