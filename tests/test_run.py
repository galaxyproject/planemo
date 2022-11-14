"""The module contains a class to test the ``cwl_run`` command."""
import os

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
                tool_path,
                job_path,
            ]
            self._check_exit_code(test_cmd)
            assert os.path.exists(os.path.join(f, "tool_test_output.html"))
            assert os.path.exists(os.path.join(f, "tool_test_output.json"))

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
                tool_path,
                job_path,
            ]
            self._check_exit_code(test_cmd)
            assert os.path.exists(os.path.join(f, "tool_test_output.html"))
            assert os.path.exists(os.path.join(f, "tool_test_output.json"))

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
