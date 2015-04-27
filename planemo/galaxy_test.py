""" Utilities for reasoning about Galaxy test results.
"""
from __future__ import absolute_import

from collections import namedtuple
import json
import xml.etree.ElementTree as ET


class GalaxyTestResults(object):
    """ Class that combine the test-centric xunit output
    with the Galaxy centric structured data output - and
    abstracts away the difference (someday).
    """

    def __init__(
        self,
        output_json_path,
        output_xml_path,
        output_html_path,
        exit_code,
    ):
        self.output_html_path = output_html_path
        self.exit_code = exit_code
        try:
            output_json_f = open(output_json_path, "r")
            structured_data = json.load(output_json_f)
            structured_data_tests = structured_data["tests"]
        except Exception:
            # Older Galaxy's will not support this option.
            structured_data = {}
            structured_data_tests = {}
        self.structured_data = structured_data
        self.structured_data_tests = structured_data_tests

        structured_data_by_id = {}
        for test in self.structured_data_tests:
            structured_data_by_id[test["id"]] = test["data"]
        self.structured_data_by_id = structured_data_by_id

        if output_xml_path:
            self.xunit_tree = parse_xunit_report(output_xml_path)
            self.__merge_xunit()
            self.has_details = True
        else:
            self.xunit_tree = ET.fromstring("<testsuite />")
            self.has_details = False
        try:
            json.dump(self.structured_data, open(output_json_path, "w"))
        except Exception:
            pass

    @property
    def _xunit_root(self):
        return self.xunit_tree.getroot()

    def __merge_xunit(self):
        xunit_attrib = self._xunit_root.attrib
        num_tests = int(xunit_attrib.get("tests", 0))
        num_failures = int(xunit_attrib.get("failures", 0))
        num_errors = int(xunit_attrib.get("errors", 0))
        num_skips = int(xunit_attrib.get("skips", 0))
        summary = dict(
            num_tests=num_tests,
            num_failures=num_failures,
            num_errors=num_errors,
            num_skips=num_skips,
        )

        self.structured_data["summary"] = summary
        self.num_tests = num_tests
        self.num_problems = num_skips + num_errors + num_failures

        for testcase_el in self.xunit_testcase_elements:
            test = case_id(testcase_el)
            test_data = self.structured_data_by_id.get(test.id)
            if not test_data:
                continue
            problem_el = None
            for problem_type in ["skip", "failure", "error"]:
                problem_el = testcase_el.find(problem_type)
                if problem_el is not None:
                    break
            if problem_el is not None:
                status = problem_el.tag
                test_data["problem_type"] = problem_el.attrib["type"]
                test_data["problem_log"] = problem_el.text
            else:
                status = "success"
            test_data["status"] = status

    @property
    def all_tests_passed(self):
        return self.num_problems == 0

    @property
    def xunit_testcase_elements(self):
        for testcase_el in find_cases(self._xunit_root):
            yield testcase_el


def parse_xunit_report(xunit_report_path):
    return ET.parse(xunit_report_path)


def find_cases(xunit_root):
    return xunit_root.findall("testcase")


def case_id(testcase_el):
    name_raw = testcase_el.attrib["name"]
    if "TestForTool_" in name_raw:
        raw_id = name_raw
    else:
        class_name = testcase_el.attrib["classname"]
        raw_id = "{0}.{1}".format(class_name, name_raw)

    name = None
    num = None
    if "TestForTool_" in raw_id:
        tool_and_num = raw_id.split("TestForTool_")[-1]
        if ".test_tool_" in tool_and_num:
            name, num_str = tool_and_num.split(".test_tool_", 1)
            num = _parse_num(num_str)
            # Tempted to but something human friendly in here like
            # num + 1 - but then it doesn't match HTML report.
        else:
            name = tool_and_num
    else:
        name = raw_id

    return TestId(name, num, raw_id)


def _parse_num(num_str):
    try:
        num = int(num_str)
    except ValueError:
        num = None
    return num


TestId = namedtuple("TestId", ["name", "num", "id"])


@property
def _label(self):
    if self.num is not None:
        return "{0}[{1}]".format(self.name, self.num)
    else:
        return self.id


TestId.label = _label
