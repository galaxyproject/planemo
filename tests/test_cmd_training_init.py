"""Tests for the ``training_init`` command."""
import os

from .test_utils import (
    CliTestCase,
    skip_if_environ,
    TEST_DATA_DIR,
)


class CmdTrainingInitTestCase(CliTestCase):
    """Container class defining test cases for the ``training_init`` command."""

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    def test_training_init_command_by_default(self):
        """Test training_init command with only topic name."""
        with self._isolate():
            training_init_command = ["training_init", "--topic_name", "test"]
            self._check_exit_code(training_init_command, exit_code=0)

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    def test_training_init_command_topic(self):
        """Test training_init command to create new topic."""
        with self._isolate():
            # working test
            training_init_command = [
                "training_init",
                "--topic_name",
                "test",
                "--topic_title",
                "Topic title",
                "--topic_target",
                "use",
                "--topic_summary",
                "Summary",
            ]
            self._check_exit_code(training_init_command, exit_code=0)
            # failing test
            training_init_command = [
                "training_init",
                "--topic_name",
                "test",
                "--topic_title",
                "Topic title",
                "--topic_target",
                "test",
                "--topic_summary",
                "Summary",
            ]
            self._check_exit_code(training_init_command, exit_code=2)

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    def test_training_init_command_tutorial_no_topic(self):
        """Test training_init command with tutorial but no topic."""
        with self._isolate():
            # working test
            training_init_command = ["training_init", "--tutorial_name", "test"]
            self._check_exit_code(training_init_command, exit_code=2)

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    def test_training_init_command_tutorial(self):
        """Test training_init command to create new tutorial."""
        with self._isolate():
            # working test
            training_init_command = [
                "training_init",
                "--topic_name",
                "test",
                "--tutorial_name",
                "test",
                "--tutorial_title",
                "Title of the tutorial",
                "--hands_on",
                "--slides",
            ]
            self._check_exit_code(training_init_command, exit_code=0)

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    def test_training_init_command_tutorial_zenodo(self):
        """Test training_init command to create new tutorial with zenodo."""
        with self._isolate():
            datatype = os.path.join(TEST_DATA_DIR, "training_datatypes.yaml")
            # not working test
            training_init_command = [
                "training_init",
                "--topic_name",
                "test",
                "--tutorial_name",
                "test",
                "--zenodo_link",
                "https://zenodo.org/record/1321885",
            ]
            self._check_exit_code(training_init_command, exit_code=1)
            # working
            training_init_command = [
                "training_init",
                "--topic_name",
                "test",
                "--tutorial_name",
                "test",
                "--zenodo_link",
                "https://zenodo.org/record/1321885",
                "--datatypes",
                datatype,
            ]
            self._check_exit_code(training_init_command, exit_code=0)

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    def test_training_init_command_tutorial_local_wf(self):
        """Test training_init command to create new tutorial with local workflow."""
        with self._isolate():
            test_workflow = os.path.join(TEST_DATA_DIR, "test_workflow_1.ga")
            # working test
            training_init_command = [
                "training_init",
                "--topic_name",
                "test",
                "--tutorial_name",
                "test",
                "--workflow",
                test_workflow,
            ]
            self._check_exit_code(training_init_command, exit_code=0)

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    def test_training_init_command_tutorial_remote_wf(self):
        """Test training_init command to create new tutorial with workflow on running instance."""
        with self._isolate():
            # not working test
            training_init_command = [
                "training_init",
                "--topic_name",
                "test",
                "--tutorial_name",
                "test",
                "--workflow_id",
                "ID",
            ]
            self._check_exit_code(training_init_command, exit_code=1)
            # working test
            training_init_command = [
                "training_init",
                "--topic_name",
                "test",
                "--tutorial_name",
                "test",
                "--workflow_id",
                "ID",
                "--galaxy_url",
                "https://usegalaxy.eu/",
                "--galaxy_api_key",
                "API",
            ]
            self._check_exit_code(training_init_command, exit_code=0)
