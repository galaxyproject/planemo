"""Tests for the ``dockstore_init`` command."""

import os
from typing import Optional

import yaml

from .test_utils import CliTestCase


class CmdDockstoreInitTestCase(CliTestCase):
    def test_plain_init(self):
        self.run_dockstore_init()

    def test_init_publish_true(self):
        self.run_dockstore_init(True)

    def test_init_publish_false(self):
        self.run_dockstore_init(False)

    def test_init_with_creator(self):
        self.run_dockstore_init_with_creator()

    def test_multiple_workflow_names_strip_complete_suffix(self):
        workflow_contents = {
            "alpha.gxwf.yml": {"class": "GalaxyWorkflow"},
            "beta.gxwf.yaml": {"class": "GalaxyWorkflow"},
            "my.gadget.ga": {"a_galaxy_workflow": "true"},
            "generic.yml": {"class": "GalaxyWorkflow"},
            "generic-yaml.yaml": {"class": "GalaxyWorkflow"},
        }
        expected_names = {
            f"/{path}": expected_name
            for path, expected_name in (
                ("alpha.gxwf.yml", "alpha"),
                ("beta.gxwf.yaml", "beta"),
                ("my.gadget.ga", "my.gadget"),
                ("generic.yml", "generic.yml"),
                ("generic-yaml.yaml", "generic-yaml.yaml"),
            )
        }

        with self._isolate():
            for path, contents in workflow_contents.items():
                with open(path, "w") as workflow_file:
                    yaml.safe_dump(contents, workflow_file)

            self._check_exit_code(["dockstore_init"])

            with open(".dockstore.yml") as dockstore_file:
                dockstore_config = yaml.safe_load(dockstore_file)
            actual_names = {
                workflow["primaryDescriptorPath"]: workflow["name"] for workflow in dockstore_config["workflows"]
            }
            assert actual_names == expected_names

    def run_dockstore_init(self, publish: Optional[bool] = None):
        with self._isolate_with_test_data("wf_repos/from_format2/0_basic_native") as f:
            init_cmd = ["dockstore_init"]
            expect_published = True
            if publish:
                init_cmd.append("--publish")
            elif publish is False:
                expect_published = False
                init_cmd.append("--no_publish")
            self._check_exit_code(init_cmd)
            assert os.path.exists(".dockstore.yml")
            with open(".dockstore.yml") as fh:
                dockstore_config = yaml.safe_load(fh)
            assert str(dockstore_config["version"]) == "1.2"
            assert "workflows" in dockstore_config
            assert len(dockstore_config["workflows"]) == 1
            assert dockstore_config["workflows"][0]["publish"] == expect_published
            workflow_lint_cmd = ["workflow_lint", "--fail_level", "error", f]
            self._check_exit_code(workflow_lint_cmd)

    def run_dockstore_init_with_creator(self):
        with self._isolate_workflow("wf_repos/autoupdate_tests/workflow_with_unexisting_tool.ga") as f:
            init_cmd = ["dockstore_init"]
            self._check_exit_code(init_cmd)
            assert os.path.exists(".dockstore.yml")
            with open(".dockstore.yml") as fh:
                dockstore_config = yaml.safe_load(fh)
            assert str(dockstore_config["version"]) == "1.2"
            assert "workflows" in dockstore_config
            assert len(dockstore_config["workflows"]) == 1
            workflow = dockstore_config["workflows"][0]
            assert "authors" in workflow
            assert workflow["authors"][0]["orcid"] == "0000-0002-1964-4960"
            assert workflow["authors"][1]["orcid"] == "0000-0002-2799-424X"
            workflow_lint_cmd = ["workflow_lint", "--skip", "tool_ids", "--fail_level", "error", f]
            self._check_exit_code(workflow_lint_cmd)
