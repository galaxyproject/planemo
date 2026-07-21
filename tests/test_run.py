"""The module contains a class to test the ``cwl_run`` command."""

import json
import os
from uuid import uuid4

import pytest

from planemo.output_models import PlanemoRunOutputs
from planemo.test.models import PlanemoTestReport
from .test_utils import (
    CliTestCase,
    CWL_DRAFT3_DIR,
    mark,
    PROJECT_TEMPLATES_DIR,
    skip_if_environ,
    target_galaxy_branch,
    TEST_DATA_DIR,
)


def _cwl_file(name):
    return os.path.normpath(os.path.join(CWL_DRAFT3_DIR, name))


# TODO: Improve these tests so they actually check something instead
# of just arbitrarily exercising the code.
class RunTestCase(CliTestCase):
    @skip_if_environ("PLANEMO_SKIP_CWLTOOL_TESTS")
    def test_run_cat_cwltool(self):
        with self._isolate() as f:
            tool_path = _cwl_file("cat1-tool.cwl")
            job_path = _cwl_file("cat-job.json")
            test_cmd = [
                "run",
                "--engine",
                "cwltool",
                "--no_container",
                "--output_json",
                "run_outputs.json",
                tool_path,
                job_path,
            ]
            self._check_exit_code(test_cmd)
            assert os.path.exists(os.path.join(f, "tool_test_output.html"))
            assert os.path.exists(os.path.join(f, "tool_test_output.json"))
            with open(os.path.join(f, "tool_test_output.json")) as test_report:
                PlanemoTestReport.model_validate(json.load(test_report))
            with open(os.path.join(f, "run_outputs.json")) as run_outputs:
                PlanemoRunOutputs.model_validate(json.load(run_outputs))

    @skip_if_environ("PLANEMO_SKIP_CWLTOOL_TESTS")
    def test_run_cat_cwltool_more_options(self):
        with self._isolate() as f:
            tool_path = _cwl_file("cat1-tool.cwl")
            job_path = _cwl_file("cat-job.json")
            test_cmd = [
                "--verbose",
                "run",
                "--engine",
                "cwltool",
                "--no_container",
                tool_path,
                job_path,
            ]
            self._check_exit_code(test_cmd)
            assert os.path.exists(os.path.join(f, "tool_test_output.html"))
            assert os.path.exists(os.path.join(f, "tool_test_output.json"))

    @skip_if_environ("PLANEMO_SKIP_CWLTOOL_TESTS")
    def test_run_cwltool_use_cache(self):
        with self._isolate() as f:
            cache_directory = os.path.join(f, "cwltool_cache")
            test_cmd = [
                "run",
                "--engine",
                "cwltool",
                "--no_container",
                "--cwltool_cache_directory",
                cache_directory,
                _cwl_file("cat1-tool.cwl"),
                _cwl_file("cat-job.json"),
            ]
            self._check_exit_code(test_cmd)
            assert os.listdir(cache_directory)
            result = self._check_exit_code(test_cmd)
            assert "Using cached output in" in result.output

    @skip_if_environ("PLANEMO_SKIP_CWLTOOL_TESTS")
    def test_run_cwltool_no_use_cache(self):
        with self._isolate() as f:
            cache_directory = os.path.join(f, "cwltool_cache")
            test_cmd = [
                "run",
                "--engine",
                "cwltool",
                "--no_container",
                "--no_use_cache",
                "--cwltool_cache_directory",
                cache_directory,
                _cwl_file("cat1-tool.cwl"),
                _cwl_file("cat-job.json"),
            ]
            self._check_exit_code(test_cmd)
            assert not os.path.exists(cache_directory)

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    @mark.tests_galaxy_branch
    def test_run_gxtool_randomlines(self):
        with self._isolate() as f:
            tool_path = os.path.join(PROJECT_TEMPLATES_DIR, "demo", "randomlines.xml")
            job_path = os.path.join(TEST_DATA_DIR, "randomlines_job_1.json")
            test_cmd = [
                "--verbose",
                "run",
                "--no_dependency_resolution",
                "--galaxy_branch",
                target_galaxy_branch(),
                "--test_data",
                TEST_DATA_DIR,
                tool_path,
                job_path,
            ]
            self._check_exit_code(test_cmd)
            assert os.path.exists(os.path.join(f, "tool_test_output.html"))
            assert os.path.exists(os.path.join(f, "tool_test_output.json"))

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    @mark.tests_galaxy_branch
    def test_run_gxtool_use_cache(self):
        # Galaxy only reuses a job if it finds an equivalent one belonging to the same
        # user - so this needs a profile to keep the database and files around between
        # runs. job_properties takes no data inputs, if it did each run would upload a
        # fresh copy of the input and Galaxy would never consider the jobs equivalent.
        tool_path = os.path.join(TEST_DATA_DIR, "tools", "functional_test_tools", "job_properties.xml")
        job_path = os.path.join(TEST_DATA_DIR, "job_properties_job.json")
        # unique so an interrupted run does not leave a profile behind that blocks the next one
        profile_name = f"planemo_test_use_cache_{uuid4().hex[:8]}"

        with self._isolate() as f:
            self._check_exit_code(["profile_create", profile_name, "--database_type", "sqlite"])
            try:

                def run(label, *extra_args):
                    report_path = os.path.join(f, f"{label}_report.json")
                    self._check_exit_code(
                        [
                            "run",
                            "--no_dependency_resolution",
                            "--galaxy_branch",
                            target_galaxy_branch(),
                            "--profile",
                            profile_name,
                            "--test_output_json",
                            report_path,
                            *extra_args,
                            tool_path,
                            job_path,
                        ]
                    )
                    with open(report_path) as report_f:
                        return json.load(report_f)["tests"][0]["data"]["job"]

                first_job = run("first")
                assert first_job["copied_from_job_id"] is None, first_job
                cached_job = run("cached")
                assert cached_job["copied_from_job_id"] is not None, cached_job
                uncached_job = run("uncached", "--no_use_cache")
                assert uncached_job["copied_from_job_id"] is None, uncached_job
            finally:
                self._check_exit_code(["profile_delete", profile_name])

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    @skip_if_environ("PLANEMO_SKIP_CWLTOOL_TESTS")
    def test_run_cat(self):
        with self._isolate() as f:
            tool_path = _cwl_file("cat1-tool.cwl")
            job_path = _cwl_file("cat-job.json")
            test_cmd = [
                "run",
                "--no_dependency_resolution",
                tool_path,
                job_path,
            ]
            self._check_exit_code(test_cmd)
            assert os.path.exists(os.path.join(f, "tool_test_output.html"))
            assert os.path.exists(os.path.join(f, "tool_test_output.json"))

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    @skip_if_environ("PLANEMO_SKIP_CWLTOOL_TESTS")
    @skip_if_environ("PLANEMO_SKIP_GALAXY_CWL_TESTS")
    def test_run_output_directory(self):
        with self._isolate() as f:
            tool_path = _cwl_file("wc-tool.cwl")
            job_path = _cwl_file("wc-job.json")
            test_cmd = [
                "--verbose",
                "run",
                "--no_dependency_resolution",
                "--output_directory",
                f,
                tool_path,
                job_path,
            ]
            self._check_exit_code(test_cmd)
            assert os.path.exists(os.path.join(f, "tool_test_output.html"))
            assert os.path.exists(os.path.join(f, "tool_test_output.json"))
            output_path = os.path.join(f, "output")
            assert os.path.exists(output_path)
            with open(output_path) as fh:
                assert fh.read().startswith("  16  198 1111")

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    @skip_if_environ("PLANEMO_SKIP_CWLTOOL_TESTS")
    @skip_if_environ("PLANEMO_SKIP_GALAXY_CWL_TESTS")
    def test_run_download_output(self):
        with self._isolate() as f:
            cat = os.path.join(PROJECT_TEMPLATES_DIR, "demo", "cat.xml")
            wf = os.path.join(TEST_DATA_DIR, "wf2.ga")
            job_path = os.path.join(TEST_DATA_DIR, "wf2-job.yml")
            output_path = os.path.join(f, "output")
            if not os.path.exists(output_path):
                os.makedirs(output_path)
            test_cmd = [
                "--verbose",
                "run",
                "--no_dependency_resolution",
                "--extra_tools",
                cat,
                "--download_outputs",
                "--output_directory",
                output_path,
                wf,
                job_path,
            ]
            self._check_exit_code(test_cmd)
            output_files = os.listdir(output_path)
            print(f"Files in {output_path}:", output_files)
            assert len(output_files) == 1
            assert output_files[0].endswith(".txt")

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    def test_run_export_invocation(self):
        with self._isolate() as f:
            cat = os.path.join(PROJECT_TEMPLATES_DIR, "demo", "cat.xml")
            wf = os.path.join(TEST_DATA_DIR, "wf2.ga")
            job_path = os.path.join(TEST_DATA_DIR, "wf2-job.yml")
            export_path = os.path.join(f, "invocation_export.rocrate.zip")
            test_cmd = [
                "--verbose",
                "run",
                "--no_dependency_resolution",
                "--extra_tools",
                cat,
                "--export_invocation",
                export_path,
                wf,
                job_path,
            ]
            self._check_exit_code(test_cmd)

            # Check that the export file was created
            assert os.path.exists(export_path)
            assert os.path.getsize(export_path) > 0

            # Verify it's a valid zip file by checking it can be opened
            import zipfile

            with zipfile.ZipFile(export_path, "r") as zip_ref:
                # Should contain some files for a valid RO-Crate
                assert len(zip_ref.namelist()) > 0

    @pytest.mark.skipif(
        target_galaxy_branch() == "release_22.05",
        reason="Skipping test on Galaxy 22.05, TRS import not supported.",
    )
    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    @mark.tests_galaxy_branch
    def test_run_trs_id(self):
        """Test importing and running a workflow using a TRS ID from GitHub."""
        with self._isolate() as f:
            # Use a TRS ID format: #workflow/github.com/org/repo/workflow_name[/version]
            # Testing with a simple workflow that uses the cat tool
            # This workflow exists on Dockstore and has a "master" version
            trs_id = "#workflow/github.com/jmchilton/galaxy-workflow-dockstore-example-1/mycoolworkflow"
            cat = os.path.join(PROJECT_TEMPLATES_DIR, "demo", "cat.xml")
            wf_job = os.path.join(TEST_DATA_DIR, "wf3-job.yml")

            test_cmd = [
                "--verbose",
                "run",
                "--no_dependency_resolution",
                "--extra_tools",
                cat,
                "--galaxy_branch",
                target_galaxy_branch(),
                "--test_data",
                TEST_DATA_DIR,
                trs_id,
                wf_job,
            ]
            self._check_exit_code(test_cmd)
            assert os.path.exists(os.path.join(f, "tool_test_output.html"))
            assert os.path.exists(os.path.join(f, "tool_test_output.json"))
