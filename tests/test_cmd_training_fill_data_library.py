"""Tests for the ``training_fill_data_library`` command."""
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


class CmdTrainingFillDataLibraryTestCase(CliTestCase):
    """Container class defining test cases for the ``training_fill_data_library`` command."""

    def test_training_fill_data_library_command_empty(self):
        with self._isolate():
            training_fill_data_library_command = [
                "training_fill_data_library"
            ]
            self._check_exit_code(training_fill_data_library_command, exit_code=2)

    def test_training_fill_data_library_command_topic(self):
        with self._isolate():
            training_fill_data_library_command = [
                "training_fill_data_library",
                "--topic_name", "test"
            ]
            self._check_exit_code(training_fill_data_library_command, exit_code=2)

    def test_training_fill_data_library_command_tutorial_topic(self):
        with self._isolate():
            topic_n = "test"
            tuto_n = "test"
            datatype = os.path.join(TEST_DATA_DIR, "training_datatypes.yaml")
            # not working
            create_tutorial_dir(topic_n, tuto_n, "training_metadata_wo_zenodo.yaml")
            training_fill_data_library_command = [
                "training_fill_data_library",
                "--topic_name", topic_n,
                "--tutorial_name", tuto_n,
                "--datatypes", datatype
            ]
            remove_topics()
            # working
            self._check_exit_code(training_fill_data_library_command, exit_code=-1)
            create_tutorial_dir(topic_n, tuto_n, "training_metadata_w_zenodo.yaml")
            training_fill_data_library_command = [
                "training_fill_data_library",
                "--topic_name", topic_n,
                "--tutorial_name", tuto_n,
                "--datatypes", datatype
            ]
            self._check_exit_code(training_fill_data_library_command, exit_code=0)

    def test_training_fill_data_library_command_tutorial_zenodo(self):
        with self._isolate():
            topic_n = "test"
            tuto_n = "test"
            create_tutorial_dir(topic_n, tuto_n, "training_metadata_wo_zenodo.yaml")
            datatype = os.path.join(TEST_DATA_DIR, "training_datatypes.yaml")
            # not working test
            training_fill_data_library_command = [
                "training_fill_data_library",
                "--topic_name", topic_n,
                "--tutorial_name", tuto_n,
                "--zenodo", "https://zenodo.org/record/1321885"
            ]
            self._check_exit_code(training_fill_data_library_command, exit_code=-1)
            # working
            training_fill_data_library_command = [
                "training_fill_data_library",
                "--topic_name", topic_n,
                "--tutorial_name", tuto_n,
                "--zenodo", "https://zenodo.org/record/1321885",
                "--datatypes", datatype
            ]
            self._check_exit_code(training_fill_data_library_command, exit_code=0)
