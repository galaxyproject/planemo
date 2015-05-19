import os

from .test_utils import (
    CliTestCase,
    skip_if_environ,
    TEST_REPOS_DIR,
)


class ShedTestTestCase(CliTestCase):

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    def test_serve(self):
        fastqc_path = os.path.join(TEST_REPOS_DIR, "fastqc")
        test_cmd = [
            "shed_test",
            "--install_galaxy",
            fastqc_path,
        ]
        self._check_exit_code(test_cmd)
