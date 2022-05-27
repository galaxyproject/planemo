"""Tests for planemo commands relating to external Galaxy instances
"""
import os

import yaml

from planemo import cli
from planemo.engine import engine_context
from planemo.runnable import for_path
from .test_utils import (
    CliTestCase,
    PROJECT_TEMPLATES_DIR,
    TEST_DATA_DIR,
)


class ExternalGalaxyCommandsTestCase(CliTestCase):
    def test_plain_init(self):
        ctx = cli.PlanemoCliContext()
        ctx.planemo_directory = "/tmp/planemo-test-workspace"
        cat_tool = os.path.join(PROJECT_TEMPLATES_DIR, "demo", "cat.xml")
        test_workflow_path = os.path.join(TEST_DATA_DIR, "wf2.ga")

        with engine_context(ctx, extra_tools=(cat_tool,)) as galaxy_engine:
            with galaxy_engine.ensure_runnables_served([for_path(test_workflow_path)]) as config:
                wfid = config.workflow_id(test_workflow_path)

                # commands to test
                profile_list_cmd = ["profile_list"]
                profile_create_cmd = [
                    "profile_create",
                    "test_ext_profile",
                    "--galaxy_url",
                    config.galaxy_url,
                    "--galaxy_user_key",
                    config.user_api_key,
                ]
                alias_create_cmd = ["create_alias", wfid, "--alias", "test_wf_alias", "--profile", "test_ext_profile"]
                alias_list_cmd = ["list_alias", "--profile", "test_ext_profile"]
                alias_delete_cmd = ["delete_alias", "--alias", "test_wf_alias", "--profile", "test_ext_profile"]
                profile_delete_cmd = ["profile_delete", "test_ext_profile"]
                run_cmd = [
                    "run",
                    "test_wf_alias",
                    os.path.join(TEST_DATA_DIR, "wf2-job.yml"),
                    "--profile",
                    "test_ext_profile",
                ]
                list_invocs_cmd = ["list_invocations", "test_wf_alias", "--profile", "test_ext_profile"]
                rerun_cmd = ["rerun", "--invocation", "invocation_id", "--profile", "test_ext_profile"]
                upload_data_cmd = [
                    "upload_data",
                    "test_wf_alias",
                    os.path.join(TEST_DATA_DIR, "wf2-job.yml"),
                    "new-job.yml",
                    "--profile",
                    "test_ext_profile",
                ]
                workflow_test_init_cmd = [
                    "workflow_test_init",
                    "invocation_id",
                    "--from_invocation",
                    "--profile",
                    "test_ext_profile",
                ]
                test_workflow_test_init_cmd = ["test", "TestWorkflow1.ga", "--profile", "test_ext_profile"]

                # test alias and profile creation
                result = self._check_exit_code(profile_list_cmd)
                assert "test_ext_profile" not in result.output
                result = self._check_exit_code(profile_create_cmd)
                assert "Profile [test_ext_profile] created" in result.output
                result = self._check_exit_code(profile_list_cmd)
                assert "test_ext_profile" in result.output
                result = self._check_exit_code(alias_create_cmd)
                assert "Alias test_wf_alias created." in result.output
                result = self._check_exit_code(alias_list_cmd)
                assert "test_wf_alias" in result.output
                assert wfid in result.output
                assert "1 aliases were found for profile test_ext_profile." in result.output

                # test upload_data command
                self._check_exit_code(upload_data_cmd)
                with open("new-job.yml") as f:
                    new_job = yaml.safe_load(f)
                assert list(new_job["WorkflowInput1"]) == list(new_job["WorkflowInput2"]) == ["class", "galaxy_id"]

                # test WF execution (from wfid) using created profile and alias
                result = self._check_exit_code(run_cmd)
                assert "Run failed" not in result.output
                result = self._check_exit_code(run_cmd + ["--no_wait"])
                assert "Run successfully executed" in result.output
                result = self._check_exit_code(list_invocs_cmd)
                assert "2 invocations found." in result.output
                assert (
                    "1 jobs ok" in result.output or '"ok": 1' in result.output
                )  # so it passes regardless if tabulate is installed or not

                # test rerun
                invocation_id = config.user_gi.workflows.get_invocations(wfid)[0]["id"]
                rerun_cmd[2] = invocation_id
                result = self._check_exit_code(rerun_cmd)
                assert "No jobs matching the specified invocation" in result.output

                # test generating test case from invocation_id
                workflow_test_init_cmd[1] = invocation_id
                self._check_exit_code(workflow_test_init_cmd)
                assert os.path.exists("TestWorkflow1.ga")
                assert os.path.exists("TestWorkflow1-tests.yml")
                self._check_exit_code(test_workflow_test_init_cmd)

                # test alias and profile deletion
                result = self._check_exit_code(alias_delete_cmd)
                assert "Alias test_wf_alias was successfully deleted from profile test_ext_profile" in result.output
                result = self._check_exit_code(alias_list_cmd)
                assert "0 aliases were found for profile test_ext_profile." in result.output
                result = self._check_exit_code(profile_delete_cmd)
                assert "Profile deleted." in result.output
                result = self._check_exit_code(profile_list_cmd)
                assert "test_ext_profile" not in result.output
