""" Tests for the galaxy_test module :) """

import os

from .test_utils import TEST_DATA_DIR

from planemo.galaxy.test import structures
from planemo.galaxy.test.actions import passed

nose_1_3_report = os.path.join(TEST_DATA_DIR, "xunit_nose_1_3.xml")
nose_0_11_report = os.path.join(TEST_DATA_DIR, "xunit_nose_0_11.xml")

xunit_report_with_failure = os.path.join(TEST_DATA_DIR, "xunit_failure.xml")


def get_test_id_new():
    _get_test_id(nose_1_3_report)


def get_test_id_old():
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
    xml_tree = structures.parse_xunit_report(xunit_report_with_failure)
    root = xml_tree.getroot()
    good_testcase_el = structures.find_cases(root)[0]
    assert passed(good_testcase_el)

    bad_testcase_el = structures.find_cases(root)[1]
    assert not passed(bad_testcase_el)
