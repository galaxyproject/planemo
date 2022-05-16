"""Tests for the ``training_generate_from_wf`` command."""
import os
import shutil

from .test_utils import (
    CliTestCase,
    skip_if_environ,
    TEST_DATA_DIR,
)


def create_tutorial_dir(topic_n, tuto_n):
    """Create the tutorial directory structure."""
    topic_dir = os.path.join("topics", topic_n)
    tuto_dir = os.path.join(topic_dir, "tutorials", tuto_n)
    metadata_path = os.path.join(topic_dir, "metadata.yaml")
    if not os.path.isdir(tuto_dir):
        os.makedirs(tuto_dir)
    if not os.path.exists(metadata_path):
        metadata = os.path.join(TEST_DATA_DIR, "training_metadata.yaml")
        shutil.copy(metadata, metadata_path)
    shutil.copy(os.path.join(TEST_DATA_DIR, "training_tutorial.md"), os.path.join(tuto_dir, "tutorial.md"))


class CmdTrainingGenerateFromWfTestCase(CliTestCase):
    """Container class defining test cases for the ``training_generate_from_wf`` command."""

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    def test_training_generate_from_wf_command_empty(self):
        """Test training_generate_from_wf command with no arguments."""
        with self._isolate():
            training_fill_data_library_command = ["training_generate_from_wf"]
            self._check_exit_code(training_fill_data_library_command, exit_code=2)

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    def test_training_generate_from_wf_command_topic(self):
        """Test training_generate_from_wf command with only topic name."""
        with self._isolate():
            training_fill_data_library_command = ["training_generate_from_wf", "--topic_name", "test"]
            self._check_exit_code(training_fill_data_library_command, exit_code=2)

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    def test_training_generate_from_wf_command_local_wf(self):
        """Test training_generate_from_wf command with local workflow."""
        with self._isolate():
            topic_n = "test"
            tuto_n = "test"
            test_workflow = os.path.join(TEST_DATA_DIR, "test_workflow_1.ga")
            # working test
            create_tutorial_dir(topic_n, tuto_n)
            training_init_command = [
                "training_generate_from_wf",
                "--topic_name",
                topic_n,
                "--tutorial_name",
                tuto_n,
                "--workflow",
                test_workflow,
            ]
            self._check_exit_code(training_init_command, exit_code=0)
            shutil.rmtree("topics")

    @skip_if_environ("PLANEMO_SKIP_GALAXY_TESTS")
    def test_training_generate_from_wf_command_remote_wf(self):
        """Test training_generate_from_wf command with workflow on running instance."""
        with self._isolate():
            topic_n = "test"
            tuto_n = "test"
            # not working test
            training_init_command = [
                "training_generate_from_wf",
                "--topic_name",
                topic_n,
                "--tutorial_name",
                tuto_n,
                "--workflow_id",
                "ID",
            ]
            self._check_exit_code(training_init_command, exit_code=self.non_zero_exit_code)
            # not working test
            create_tutorial_dir(topic_n, tuto_n)
            training_init_command = [
                "training_generate_from_wf",
                "--topic_name",
                topic_n,
                "--tutorial_name",
                tuto_n,
                "--workflow_id",
                "ID",
                "--galaxy_url",
                "https://usegalaxy.eu/",
                "--galaxy_api_key",
                "API",
            ]
            self._check_exit_code(training_init_command, exit_code=self.non_zero_exit_code)
            shutil.rmtree("topics")
