import os

from .test_utils import (
    CliTestCase,
    skip_if_environ,
    PROJECT_TEMPLATES_DIR,
    skip_unless_python_2_7,
)

CWL_DRAFT2_DIR = os.path.join(PROJECT_TEMPLATES_DIR, "cwl_draft2_spec")


def _cwl_file(name):
    return os.path.normpath(os.path.join(CWL_DRAFT2_DIR, name))


# TODO: Improve these tests so they actually check something instead
# of just arbitrarily exercising the code.
class CwlRunTestCase(CliTestCase):

    @skip_unless_python_2_7()
    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    def test_run_cat(self):
        with self._isolate():
            tool_path = _cwl_file("cat1-tool.cwl")
            job_path = _cwl_file("cat-job.json")
            test_cmd = [
                "cwl_run",
                tool_path,
                job_path,
            ]
            self._check_exit_code(test_cmd)

    @skip_unless_python_2_7()
    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    def test_run_cat_conformance(self):
        with self._isolate():
            tool_path = _cwl_file("cat1-tool.cwl")
            job_path = _cwl_file("cat-job.json")
            test_cmd = [
                "cwl_run",
                "--conformance-test",
                tool_path,
                job_path,
            ]
            self._check_exit_code(test_cmd)
