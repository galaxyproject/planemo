"""Tests for the ``training_fill_data_library`` command."""
import os
import shutil

from .test_cmd_training_generate_from_wf import create_tutorial_dir
from .test_utils import (
    CliTestCase,
    TEST_DATA_DIR,
)


class CmdTrainingFillDataLibraryTestCase(CliTestCase):
    """Container class defining test cases for the ``training_fill_data_library`` command."""

    def test_training_fill_data_library_command_empty(self):
        """Test training_fill_data_library command with no arguments."""
        with self._isolate():
            training_fill_data_library_command = ["training_fill_data_library"]
            self._check_exit_code(training_fill_data_library_command, exit_code=2)

    def test_training_fill_data_library_command_topic(self):
        """Test training_fill_data_library command with only topic name."""
        with self._isolate():
            training_fill_data_library_command = ["training_fill_data_library", "--topic_name", "test"]
            self._check_exit_code(training_fill_data_library_command, exit_code=2)

    def test_training_fill_data_library_command_tutorial_topic(self):
        """Test training_fill_data_library command with tutorial name."""
        with self._isolate():
            topic_n = "test"
            tuto_n = "test"
            datatype = os.path.join(TEST_DATA_DIR, "training_datatypes.yaml")
            # not working
            create_tutorial_dir(topic_n, tuto_n)
            training_fill_data_library_command = [
                "training_fill_data_library",
                "--topic_name",
                topic_n,
                "--tutorial_name",
                tuto_n,
                "--datatypes",
                datatype,
            ]
            shutil.rmtree("topics")
            self._check_exit_code(training_fill_data_library_command, exit_code=1)
            # working
            create_tutorial_dir(topic_n, tuto_n)
            training_fill_data_library_command = [
                "training_fill_data_library",
                "--topic_name",
                topic_n,
                "--tutorial_name",
                tuto_n,
                "--datatypes",
                datatype,
            ]
            self._check_exit_code(training_fill_data_library_command, exit_code=0)

    def test_training_fill_data_library_command_tutorial_zenodo(self):
        """Test training_fill_data_library command with zenodo link."""
        with self._isolate():
            topic_n = "test"
            tuto_n = "test"
            create_tutorial_dir(topic_n, tuto_n)
            datatype = os.path.join(TEST_DATA_DIR, "training_datatypes.yaml")
            # not working test
            training_fill_data_library_command = [
                "training_fill_data_library",
                "--topic_name",
                topic_n,
                "--tutorial_name",
                tuto_n,
                "--zenodo_link",
                "https://zenodo.org/record/1321885",
            ]
            self._check_exit_code(training_fill_data_library_command, exit_code=1)
            # working
            training_fill_data_library_command = [
                "training_fill_data_library",
                "--topic_name",
                topic_n,
                "--tutorial_name",
                tuto_n,
                "--zenodo_link",
                "https://zenodo.org/record/1321885",
                "--datatypes",
                datatype,
            ]
            self._check_exit_code(training_fill_data_library_command, exit_code=0)
