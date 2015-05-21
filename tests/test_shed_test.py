import os

from .test_utils import (
    CliTestCase,
    skip_if_environ,
    TEST_REPOS_DIR,
)

SHED_TARGET = "testtoolshed"


class ShedTestTestCase(CliTestCase):

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    def test_serve(self):
        fastqc_path = os.path.join(TEST_REPOS_DIR, "fastqc")
        test_cmd = [
            "shed_test",
            "--shed_target", SHED_TARGET,
            "--install_galaxy",
            fastqc_path,
        ]
        self._check_exit_code(test_cmd)
