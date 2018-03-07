from .test_utils import (
    CliTestCase,
    skip_if_environ,
)


class InitAndTestTestCase(CliTestCase):

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    def test_init_and_test_master(self):
        self.__run_commands()

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    def test_init_and_test_dev(self):
        self.__run_commands(test_args=["--galaxy_branch", "dev"])

    def __run_commands(self, test_args=[]):
        with self._isolate():
            init_cmd = ["project_init", "--template", "demo", "basic"]
            self._check_exit_code(init_cmd)
            test_cmd = ["test", "--install_galaxy"] + test_args
            test_cmd += ["basic/cat.xml"]
            self._check_exit_code(test_cmd)
