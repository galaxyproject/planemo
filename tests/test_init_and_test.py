from .test_utils import (
    CliTestCase,
    skip_if_environ,
)


class InitAndTestTestCase(CliTestCase):

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    def test_init_and_test(self):
        with self._isolate():
            init_cmd = ["project_init", "--template", "demo", "basic"]
            self._check_exit_code(init_cmd)
            test_cmd = [
                "test",
                "--install_galaxy",
                "basic/cat.xml"
            ]
            self._check_exit_code(test_cmd)
