"""Module contains :class:`CmdTestTestCase` - integration tests for the ``test`` command."""
import os

from .test_utils import (
    CliTestCase,
    mark,
    PROJECT_TEMPLATES_DIR,
    skip_if_environ,
    target_galaxy_branch,
    TEST_REPOS_DIR,
    TEST_TOOLS_DIR,
)


class CmdTestCondaTestCase(CliTestCase):
    """Integration tests for the ``test`` command."""

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    @mark.tests_galaxy_branch
    def test_conda_dependencies_by_default(self):
        with self._isolate():
            bwa_test = os.path.join(PROJECT_TEMPLATES_DIR, "conda_testing", "bwa.xml")
            test_command = [
                "--verbose",
                "test",
                "--galaxy_branch",
                target_galaxy_branch(),
                bwa_test,
            ]
            self._check_exit_code(test_command, exit_code=0)

    @skip_if_environ("PLANEMO_SKIP_REDUNDANT_TESTS")  # same code path as test above.
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
    @mark.tests_galaxy_branch
    def test_conda_dependencies_version(self):
        """Test tool with wrong version and ensure it fails."""
        with self._isolate():
            # Try a failing test to ensure the primary test above isn't just passing spuriously.
            bwa_test = os.path.join(TEST_TOOLS_DIR, "bwa_wrong_version.xml")
            test_command = [
                "--verbose",
                "test",
                "--galaxy_branch",
                target_galaxy_branch(),
                "--conda_dependency_resolution",
                "--conda_auto_install",
                "--conda_auto_init",
                bwa_test,
            ]
            self._check_exit_code(test_command, exit_code=1)

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    @mark.tests_galaxy_branch
    def test_local_conda_dependencies_version(self):
        """Test a tool that requires local package builds."""
        with self._isolate() as isolated_dir:
            conda_prefix = os.path.join(isolated_dir, "miniconda3")
            fleeqtk_recipe = os.path.join(PROJECT_TEMPLATES_DIR, "conda_answers", "exercise_2", "fleeqtk")
            build_command = [
                "conda_build",
                "--conda_prefix",
                conda_prefix,
                fleeqtk_recipe,
            ]
            self._check_exit_code(build_command)
            fleeqtk_tool = os.path.join(TEST_REPOS_DIR, "conda_exercises_fleeqtk", "fleeqtk_seq.xml")
            conda_install_command = [
                "conda_install",
                "--conda_prefix",
                conda_prefix,
                "--conda_use_local",
                fleeqtk_tool,
            ]
            self._check_exit_code(conda_install_command)
            test_command = [
                "test",
                "--galaxy_branch",
                target_galaxy_branch(),
                "--conda_prefix",
                conda_prefix,
                fleeqtk_tool,
            ]
            self._check_exit_code(test_command)
