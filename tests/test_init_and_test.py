from .test_utils import (
    CliTestCase,
    mark,
    skip_if_environ,
    target_galaxy_branch,
)


class InitAndTestTestCase(CliTestCase):
    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    @mark.tests_galaxy_branch
    def test_init_and_test(self):
        with self._isolate():
            init_cmd = ["project_init", "--template", "demo", "basic"]
            self._check_exit_code(init_cmd)
            test_cmd = [
                "test",
                "--no_dependency_resolution",
                "--galaxy_branch",
                target_galaxy_branch(),
                "basic/cat.xml",
            ]
            self._check_exit_code(test_cmd)
