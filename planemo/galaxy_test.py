""" Utilities for reasoning about Galaxy test results.
"""
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
            self.xunit_tree = ET.parse(output_xml_path)
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
            id = testcase_el.get("name")
            test_data = self.structured_data_by_id.get(id)
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
        for testcase_el in self._xunit_root.findall("testcase"):
            yield testcase_el
