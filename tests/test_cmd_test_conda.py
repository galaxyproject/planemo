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
    def test_conda_dependencies_explicit_resolution(self):
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
    def test_conda_dependencies_implicit_resolution(self):
        # Shouldn't need to explicitly declare resolution of these things anymore.
        with self._isolate():
            bwa_test = os.path.join(PROJECT_TEMPLATES_DIR, "conda_testing", "bwa.xml")
            test_command = [
                "--verbose",
                "test",
                "--galaxy_branch", "dev",
                bwa_test,
            ]
            self._check_exit_code(test_command, exit_code=0)

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    def test_conda_dependencies_galaxy_16_10(self):
        # Conda resolution was changed quite a bit after this - just ensure things are still working.
        with self._isolate():
            bwa_test = os.path.join(PROJECT_TEMPLATES_DIR, "conda_testing", "bwa.xml")
            test_command = [
                "--verbose",
                "test",
                "--galaxy_branch", "release_16.10",
                bwa_test,
            ]
            self._check_exit_code(test_command, exit_code=0)

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    def test_conda_dependencies_galaxy_16_07(self):
        # Conda resolution was changed quite a bit after this - just ensure things are still working.
        with self._isolate():
            bwa_test = os.path.join(PROJECT_TEMPLATES_DIR, "conda_testing", "bwa.xml")
            test_command = [
                "--verbose",
                "test",
                "--galaxy_branch", "release_16.07",
                bwa_test,
            ]
            self._check_exit_code(test_command, exit_code=0)

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    def test_conda_dependencies_verify_branch_testing_properly(self):
        # This tries to test Conda dependency resolution with a pre-Conda Galaxy and
        # expects a failure. This verifies the other test cases are in fact testing
        # different versions of Galaxy as expected.
        with self._isolate():
            bwa_test = os.path.join(PROJECT_TEMPLATES_DIR, "conda_testing", "bwa.xml")
            test_command = [
                "--verbose",
                "test",
                "--galaxy_branch", "release_15.10",
                bwa_test,
            ]
            self._check_exit_code(test_command, exit_code=self.non_zero_exit_code)

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
