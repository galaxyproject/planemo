"""Module contains :class:`CmdTestTestCase` - integration tests for the ``test`` command."""
from .test_utils import (
    CliTestCase,
    run_verbosely,
    skip_if_environ,
)


class CmdShedTestTestCase(CliTestCase):
    """Integration tests for the ``shed_test`` command."""

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    def test_simple_shed_test(self):
        """Test --galaxy_startup_timeout."""
        with self._isolate():
            test_command = self._test_command(
                "--no_dependency_resolution",
                "--owner",
                "iuc",
                "--name",
                "collection_element_identifiers",
                "--shed_target",
                "toolshed",
            )
            self._check_exit_code(test_command, exit_code=0)

    def _test_command(self, *args):
        test_cmd = ["--verbose"] if run_verbosely() else []
        test_cmd.append("shed_test")
        test_cmd.extend(args)
        return test_cmd
