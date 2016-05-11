"""Tests for the `planeo.galaxy.test` module."""

import os
import shutil

from .test_utils import (
    test_context,
    TempDirectoryTestCase,
    TEST_DATA_DIR,
)

from planemo.galaxy.test import structures
from planemo.galaxy.test.actions import passed
from planemo.galaxy.test.actions import run_in_config


nose_1_3_report = os.path.join(TEST_DATA_DIR, "xunit_nose_1_3.xml")
nose_0_11_report = os.path.join(TEST_DATA_DIR, "xunit_nose_0_11.xml")

xunit_report_with_failure = os.path.join(TEST_DATA_DIR, "xunit_failure.xml")


class RunInConfigTestCase(TempDirectoryTestCase):
    """Test cases for ``run_in_config``."""

    def setUp(self):
        """Setup mock keywords, context, and Galaxy config for tests."""
        super(RunInConfigTestCase, self).setUp()
        td = self.temp_directory
        self.ctx = test_context()
        self.config = _MockConfig(td)
        self.kwds = {
            "test_output": os.path.join(td, "tests.html"),
            "test_output_json": os.path.join(td, "tests.json"),
            "test_output_xunit": os.path.join(td, "tests.xml"),
        }

    def test_normal_execution(self):
        """Test normal operation of run_in_config."""
        def mock_galaxy_run(ctx_, command, env, action):
            assert ctx_ is self.ctx
            assert env["test_key"] == "test_value"
            self._copy_good_artifacts(["xml", "html", "json"])
            return 0

        assert self._do_run(mock_galaxy_run) == 0

    def test_failed_to_produce_xunit(self):
        """Test an exception is thrown if not XUnit report is produced."""
        def mock_galaxy_run(ctx_, command, env, action):
            self._copy_good_artifacts(["json", "html"])
            return 0

        with self.assertRaises(Exception):
            self._do_run(mock_galaxy_run)

    def test_failed_to_produce_json(self):
        """Test an exception is thrown if not XUnit report is produced."""
        def mock_galaxy_run(ctx_, command, env, action):
            self._copy_good_artifacts(["xml", "html"])
            return 0

        with self.assertRaises(Exception):
            self._do_run(mock_galaxy_run)

    def _copy_good_artifacts(self, extensions):
        for extension in extensions:
            source = os.path.join(TEST_DATA_DIR, "tt_success.%s" % extension)
            destination = os.path.join(self.temp_directory, "tests.%s" % extension)
            shutil.copy(source, destination)

    def _do_run(self, mock_run_function):
        return run_in_config(self.ctx, self.config, run=mock_run_function, **self.kwds)


def get_test_id_new():
    """Test ID parsing on newer nose dependency."""
    _get_test_id(nose_1_3_report)


def get_test_id_old():
    """Test ID parsing on older nose dependency."""
    _get_test_id(nose_0_11_report)


def _get_test_id(path):
    xml_tree = structures.parse_xunit_report(path)
    root = xml_tree.getroot()
    first_testcase = structures.find_cases(root)[0]
    test_id = structures.case_id(first_testcase)
    assert test_id.label == "cat[0]"
    expected_id = "functional.test_toolbox.TestForTool_cat.test_tool_000000"
    assert test_id.id == expected_id
    assert test_id.num == 0


def test_passed():
    """Test :func:`passed`."""
    xml_tree = structures.parse_xunit_report(xunit_report_with_failure)
    root = xml_tree.getroot()
    good_testcase_el = structures.find_cases(root)[0]
    assert passed(good_testcase_el)

    bad_testcase_el = structures.find_cases(root)[1]
    assert not passed(bad_testcase_el)


class _MockConfig(object):

    def __init__(self, temp_directory):
        self.config_directory = temp_directory
        self.env = {"test_key": "test_value"}
        self.galaxy_root = os.path.join(self.config_directory, "galaxy-dev")
