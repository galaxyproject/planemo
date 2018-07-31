"""Tests for the ``training_init`` command."""
import os

from .test_utils import (
    CliTestCase,
    PROJECT_TEMPLATES_DIR,
    TEST_DATA_DIR
)


class CmdTrainingInitTestCase(CliTestCase):
    """Container class defining test cases for the ``training_init`` command."""

    def test_training_init_command_by_default(self):
        with self._isolate():
            training_init_command = [
                "training_init",
                "--topic_name", "test"
            ]
            self._check_exit_code(training_init_command, exit_code=-1)

    def test_training_init_command_templates(self):
        with self._isolate():
            training_template = os.path.join(PROJECT_TEMPLATES_DIR, "training")
            training_init_command = [
                "training_init",
                "--topic_name", "test",
                "--templates", training_template
            ]
            self._check_exit_code(training_init_command, exit_code=0)

    def test_training_init_command_topic(self):
        with self._isolate():
            training_template = os.path.join(PROJECT_TEMPLATES_DIR, "training")
            # working test
            training_init_command = [
                "training_init",
                "--topic_name", "test",
                "--templates", training_template,
                "--topic_title", "Topic title",
                "--topic_target", "use",
                "--topic_summary", "Summary"
            ]
            self._check_exit_code(training_init_command, exit_code=0)
            # failing test
            training_init_command = [
                "training_init",
                "--topic_name", "test",
                "--templates", training_template,
                "--topic_title", "Topic title",
                "--topic_target", "test",
                "--topic_summary", "Summary"
            ]
            self._check_exit_code(training_init_command, exit_code=2)

    def test_training_init_command_tutorial_no_topic(self):
        with self._isolate():
            training_template = os.path.join(PROJECT_TEMPLATES_DIR, "training")
            # working test
            training_init_command = [
                "training_init",
                "--tutorial_name", "test",
                "--templates", training_template,
            ]
            self._check_exit_code(training_init_command, exit_code=2)

    def test_training_init_command_tutorial(self):
        with self._isolate():
            training_template = os.path.join(PROJECT_TEMPLATES_DIR, "training")
            # working test
            training_init_command = [
                "training_init",
                "--topic_name", "test",
                "--tutorial_name", "test",
                "--templates", training_template,
                "--tutorial_title", "Title of the tutorial",
                "--hands_on",
                "--slides"
            ]
            self._check_exit_code(training_init_command, exit_code=0)

    def test_training_init_command_tutorial_zenodo(self):
        with self._isolate():
            training_template = os.path.join(PROJECT_TEMPLATES_DIR, "training")
            datatype = os.path.join(TEST_DATA_DIR, "training_datatypes.yaml")
            # not working test
            training_init_command = [
                "training_init",
                "--topic_name", "test",
                "--tutorial_name", "test",
                "--zenodo", "https://zenodo.org/record/1321885",
                "--templates", training_template
            ]
            self._check_exit_code(training_init_command, exit_code=-1)
            # working
            training_init_command = [
                "training_init",
                "--topic_name", "test",
                "--tutorial_name", "test",
                "--zenodo", "https://zenodo.org/record/1321885",
                "--datatypes", datatype,
                "--templates", training_template
            ]
            self._check_exit_code(training_init_command, exit_code=0)

    def test_training_init_command_tutorial_local_wf(self):
        with self._isolate():
            training_template = os.path.join(PROJECT_TEMPLATES_DIR, "training")
            test_workflow = os.path.join(TEST_DATA_DIR, "test_workflow_1.ga")
            # working test
            training_init_command = [
                "training_init",
                "--topic_name", "test",
                "--tutorial_name", "test",
                "--workflow", test_workflow,
                "--templates", training_template
            ]
            self._check_exit_code(training_init_command, exit_code=0)

    def test_training_init_command_tutorial_remote_wf(self):
        with self._isolate():
            training_template = os.path.join(PROJECT_TEMPLATES_DIR, "training")
            # not working test
            training_init_command = [
                "training_init",
                "--topic_name", "test",
                "--tutorial_name", "test",
                "--workflow_id", "ID",
                "--templates", training_template
            ]
            self._check_exit_code(training_init_command, exit_code=-1)
            # working test
            training_init_command = [
                "training_init",
                "--topic_name", "test",
                "--tutorial_name", "test",
                "--workflow_id", "ID",
                "--galaxy_url", "https://usegalaxy.eu/",
                "--galaxy_api_key", "API",
                "--templates", training_template
            ]
            self._check_exit_code(training_init_command, exit_code=0)
