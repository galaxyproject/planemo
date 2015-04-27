""" Tests for the galaxy_test module :) """

import os

from .test_utils import TEST_DATA_DIR

from planemo import galaxy_test

nose_1_3_report = os.path.join(TEST_DATA_DIR, "xunit_nose_1_3.xml")
nose_0_11_report = os.path.join(TEST_DATA_DIR, "xunit_nose_0_11.xml")


def get_test_id_new():
    _get_test_id(nose_1_3_report)


def get_test_id_old():
    _get_test_id(nose_0_11_report)


def _get_test_id(path):
    xml_tree = galaxy_test.parse_xunit_report(path)
    root = xml_tree.getroot()
    first_testcase = galaxy_test.find_testcases(root)[0]
    test_id = galaxy_test.test_id(first_testcase)
    assert test_id.label == "cat[0]"
    expected_id = "functional.test_toolbox.TestForTool_cat.test_tool_000000"
    assert test_id.id == expected_id
    assert test_id.num == 0
