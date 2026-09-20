import os

from .test_utils import (
    CliTestCase,
    skip,
    skip_if_environ,
    TEST_REPOS_DIR,
)

SHED_TARGET = "toolshed"


class ShedTestTestCase(CliTestCase):
    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    @skip_if_environ("PLANEMO_SKIP_SHED_TESTS")
    def test_shed_test_fastqc(self):
        self.__run_shed_test("fastqc")

    @skip  # Skip for now since toolshed has problem with op column.
    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    @skip_if_environ("PLANEMO_SKIP_SHED_TESTS")
    def test_shed_test_datamash(self):
        self.__run_shed_test("datamash")

    def __run_shed_test(self, repo_name, target=SHED_TARGET):
        repo_path = os.path.join(TEST_REPOS_DIR, repo_name)
        test_cmd = [
            "shed_test",
            "--shed_target",
            target,
            "--install_galaxy",
            repo_path,
        ]
        self._check_exit_code(test_cmd)
