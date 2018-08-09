"""Tests for the ``training_generate_from_wf`` command."""
import os
import shutil

from .test_utils import (
    CliTestCase,
    TEST_DATA_DIR
)


def create_tutorial_dir(topic_n, tuto_n, metadata_n):
    topic_dir = os.path.join("topics", topic_n)
    tuto_dir = os.path.join(topic_dir, "tutorials", tuto_n)
    metadata_path = os.path.join(topic_dir, "metadata.yaml")
    if not os.path.isdir(topic_dir):
        os.makedirs(topic_dir)
    if not os.path.isdir(tuto_dir):
        os.makedirs(tuto_dir)
    if not os.path.exists(metadata_path):
        metadata = os.path.join(TEST_DATA_DIR, metadata_n)
        shutil.copy(metadata, metadata_path)


def remove_topics():
    shutil.rmtree("topics")


class CmdTrainingGenerateFromWfTestCase(CliTestCase):
    """Container class defining test cases for the ``training_generate_from_wf`` command."""

    def test_training_generate_from_wf_command_empty(self):
        with self._isolate():
            training_fill_data_library_command = [
                "training_generate_from_wf"
            ]
            self._check_exit_code(training_fill_data_library_command, exit_code=2)

    def test_training_generate_from_wf_command_topic(self):
        with self._isolate():
            training_fill_data_library_command = [
                "training_generate_from_wf",
                "--topic_name", "test"
            ]
            self._check_exit_code(training_fill_data_library_command, exit_code=2)

    def test_training_generate_from_wf_command_local_wf(self):
        with self._isolate():
            topic_n = "test"
            tuto_n = "test"
            test_workflow = os.path.join(TEST_DATA_DIR, "test_workflow_1.ga")
            # working test
            create_tutorial_dir(topic_n, tuto_n, "training_metadata_wo_zenodo.yaml")
            training_init_command = [
                "training_generate_tuto_from_wf",
                "--topic_name", "test",
                "--tutorial_name", "test",
                "--workflow", test_workflow
            ]
            self._check_exit_code(training_init_command, exit_code=-1)
            remove_topics()

    def test_training_generate_from_wf_command_remote_wf(self):
        with self._isolate():
            topic_n = "test"
            tuto_n = "test"
            # not working test
            training_init_command = [
                "training_generate_from_wf",
                "--topic_name", "test",
                "--tutorial_name", "test",
                "--workflow_id", "ID"
            ]
            self._check_exit_code(training_init_command, exit_code=-1)
            # not working test
            create_tutorial_dir(topic_n, tuto_n, "training_metadata_wo_zenodo.yaml")
            training_init_command = [
                "training_generate_from_wf",
                "--topic_name", "test",
                "--tutorial_name", "test",
                "--workflow_id", "ID",
                "--galaxy_url", "https://usegalaxy.eu/",
                "--galaxy_api_key", "API"
            ]
            self._check_exit_code(training_init_command, exit_code=0)
            remove_topics()
