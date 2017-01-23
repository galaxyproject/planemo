"""Module contains :class:`CmdTestTestCase` - integration tests for the ``test`` command."""
import os

from .test_utils import (
    CliTestCase,
    PROJECT_TEMPLATES_DIR,
    skip_if_environ,
    TEST_TOOLS_DIR,
)


class CmdTestCondaTestCase(CliTestCase):
    """Integration tests for the ``test`` command."""

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    def test_conda_dependencies(self):
        with self._isolate():
            bwa_test = os.path.join(PROJECT_TEMPLATES_DIR, "conda_testing", "bwa.xml")
            test_command = [
                "--verbose",
                "test",
                "--conda_dependency_resolution",
                "--conda_auto_install",
                "--conda_auto_init",
                bwa_test,
            ]
            self._check_exit_code(test_command, exit_code=0)

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    def test_conda_dependencies_version(self):
        """Test testing a simple workflow with Galaxy."""
        with self._isolate():
            # Try a failing test to ensure the primary test above isn't just passing spuriously.
            bwa_test = os.path.join(TEST_TOOLS_DIR, "bwa_wrong_version.xml")
            test_command = [
                "--verbose",
                "test",
                "--conda_dependency_resolution",
                "--conda_auto_install",
                "--conda_auto_init",
                bwa_test,
            ]
            self._check_exit_code(test_command, exit_code=1)
