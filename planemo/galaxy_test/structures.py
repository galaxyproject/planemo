""" Utilities for reasoning about Galaxy test results.
"""
from __future__ import print_function
from __future__ import absolute_import

from collections import namedtuple
import json
import xml.etree.ElementTree as ET

RUN_TESTS_CMD = (
    "sh run_tests.sh --report_file %s %s %s %s"
)


class GalaxyTestCommand(object):

    def __init__(
        self,
        html_report_file,
        xunit_report_file,
        structured_report_file,
        failed=False,
        installed=False,
    ):
        self.html_report_file = html_report_file
        self.xunit_report_file = xunit_report_file
        self.structured_report_file = structured_report_file
        self.failed = failed
        self.installed = installed

    def build(self):
        xunit_report_file = self.xunit_report_file
        sd_report_file = self.structured_report_file
        html_report_file = self.html_report_file
        if xunit_report_file:
            xunit_arg = "--xunit_report_file %s" % xunit_report_file
        else:
            xunit_arg = ""
        if sd_report_file:
            sd_arg = "--structured_data_report_file %s" % sd_report_file
        else:
            sd_arg = ""
        if self.installed:
            tests = "-installed"
        else:
            tests = "functional.test_toolbox"
            if self.failed:
                sd = StructuredData(self.structured_report_file)
                failed_ids = sd.failed_ids
                tests = " ".join(failed_ids)
        return RUN_TESTS_CMD % (html_report_file, xunit_arg, sd_arg, tests)


class StructuredData(object):
    """ Abstraction around Galaxy's structured test data output.
    """

    def __init__(self, json_path):
        self.json_path = json_path
        try:
            with open(json_path, "r") as output_json_f:
                structured_data = json.load(output_json_f)
                structured_data_tests = structured_data["tests"]
        except Exception:
            print("Warning: Targetting older Galaxy which did not "
                  "produce a structured test results files.")
            structured_data = {}
            structured_data_tests = {}
        self.structured_data = structured_data
        self.structured_data_tests = structured_data_tests
        structured_data_by_id = {}
        for test in self.structured_data_tests:
            structured_data_by_id[test["id"]] = test["data"]
        self.structured_data_by_id = structured_data_by_id
        self.has_details = "summary" in structured_data
        if self.has_details:
            self._read_summary()

    def update(self):
        with open(self.json_path, "w") as out_f:
            json.dump(self.structured_data, out_f)

    def merge_xunit(self, xunit_root):
        self.has_details = True
        xunit_attrib = xunit_root.attrib
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

        for testcase_el in xunit_t_elements_from_root(xunit_root):
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

    def _read_summary(self):
        # TODO: read relevant data out of summary object.
        pass

    @property
    def failed_ids(self):
        ids = set([])
        for test_data in self.structured_data_tests:
            if test_data["data"]["status"] == "success":
                continue
            test_case = test_data["id"].replace(".test_toolbox.", ".test_toolbox:")
            ids.add(test_case)
        return ids


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
        sd = StructuredData(output_json_path)
        self.sd = sd
        self.structured_data = sd.structured_data
        self.structured_data_tests = sd.structured_data_tests
        self.structured_data_by_id = sd.structured_data_by_id

        if output_xml_path:
            self.xunit_tree = parse_xunit_report(output_xml_path)
            sd.merge_xunit(self._xunit_root)
        else:
            self.xunit_tree = ET.fromstring("<testsuite />")
        self.sd.update()

    @property
    def has_details(self):
        return self.sd.has_details

    @property
    def num_tests(self):
        return self.sd.num_tests

    @property
    def num_problems(self):
        return self.sd.num_problems

    @property
    def _xunit_root(self):
        return self.xunit_tree.getroot()

    @property
    def all_tests_passed(self):
        return self.sd.num_problems == 0

    @property
    def xunit_testcase_elements(self):
        return xunit_t_elements_from_root(self._xunit_root)


def xunit_t_elements_from_root(xunit_root):
    for testcase_el in find_cases(xunit_root):
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
