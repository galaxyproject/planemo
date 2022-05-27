""" Integration tests for shed_diff command.
"""

import os
import sys
import tempfile
from os.path import join
from xml.etree import ElementTree

from planemo import io
from planemo.xml.diff import diff
from .test_shed_upload import update_package_1
from .test_utils import (
    CliShedTestCase,
    TEST_REPOS_DIR,
)

DIFF_LINES = [
    "diff -r _workingdir_/related_file _custom_shed_/related_file",
    "< A related non-tool file (modified).",
    "> A related non-tool file.",
]


class ShedDiffTestCase(CliShedTestCase):
    def test_shed_diff(self):
        with self._isolate_repo("single_tool") as f:
            upload_command = ["shed_upload", "--force_repository_creation"]
            upload_command.extend(self._shed_args())
            self._check_exit_code(upload_command)
            io.write_file(
                join(f, "related_file"),
                "A related non-tool file (modified).\n",
            )
            self._check_diff(f, True)
            self._check_diff(f, False)

    def test_diff_doesnt_exist(self):
        with self._isolate_repo("multi_repos_nested"):
            diff_command = ["shed_diff"]
            diff_command.extend(self._shed_args(read_only=True))
            self._check_exit_code(diff_command, exit_code=2)

    def test_diff_recursive(self):
        with self._isolate_repo("multi_repos_nested") as f:
            self._shed_create(recursive=True)

            diff_command = ["shed_diff", "-r"]
            diff_command.extend(self._shed_args(read_only=True))
            self._check_exit_code(diff_command, exit_code=0)

            io.write_file(
                join(f, "cat1", "related_file"),
                "A related non-tool file (modified).\n",
            )
            self._check_exit_code(diff_command, exit_code=1)

    def test_shed_diff_raw(self):
        with self._isolate_repo("suite_auto"):
            self._shed_create()

            diff_command = ["shed_diff", "-o", "diff"]
            diff_command.append("--raw")
            diff_command.extend(self._shed_args(read_only=True))
            self._check_exit_code(diff_command, exit_code=1)

            diff_command = ["shed_diff", "-o", "diff"]
            diff_command.extend(self._shed_args(read_only=True))
            self._check_exit_code(diff_command, exit_code=0)

    def test_shed_diff_xml_no_diff(self):
        with self._isolate_repo("package_1"):
            self._shed_create()

            diff_command = ["shed_diff"]
            diff_command.extend(self._shed_args(read_only=True))
            self._check_exit_code(diff_command, exit_code=0)

    def test_shed_diff_xml_diff(self):
        with self._isolate_repo("package_1") as f:
            self._shed_create()

            update_package_1(f)

            diff_command = ["shed_diff"]
            diff_command.extend(self._shed_args(read_only=True))
            self._check_exit_code(diff_command, exit_code=1)

    def test_diff_xunit(self):
        with self._isolate_repo("multi_repos_nested") as f:
            self._shed_create(recursive=True)

            xunit_report = tempfile.NamedTemporaryFile(delete=False)
            xunit_report.flush()
            xunit_report.close()
            diff_command = ["shed_diff", "-r", "--report_xunit", xunit_report.name]
            diff_command.extend(self._shed_args(read_only=True))
            known_good_xunit_report = os.path.join(TEST_REPOS_DIR, "multi_repos_nested.xunit.xml")
            known_bad_xunit_report = os.path.join(TEST_REPOS_DIR, "multi_repos_nested.xunit-bad.xml")
            self._check_exit_code(diff_command, exit_code=0)

            with open(xunit_report.name) as fh:
                compare = fh.read()
            if diff(
                self._make_deterministic(ElementTree.parse(known_good_xunit_report).getroot()),
                self._make_deterministic(ElementTree.fromstring(compare)),
                reporter=sys.stdout.write,
            ):
                sys.stdout.write(compare)
                assert False, "XUnit report different from multi_repos_nested.xunit.xml."

            io.write_file(
                join(f, "cat1", "related_file"),
                "A related non-tool file (modified).\n",
            )
            self._check_exit_code(diff_command, exit_code=1)

            with open(xunit_report.name) as fh:
                compare = fh.read()
            if diff(
                self._make_deterministic(ElementTree.parse(known_bad_xunit_report).getroot()),
                self._make_deterministic(ElementTree.fromstring(compare)),
                reporter=sys.stdout.write,
            ):
                sys.stdout.write(compare)
                assert False, "XUnit report different from multi_repos_nested.xunit-bad.xml."

            os.unlink(xunit_report.name)

    def _check_diff(self, f, raw):
        diff_command = ["shed_diff", "-o", "diff"]
        if raw:
            diff_command.append("--raw")
        diff_command.extend(self._shed_args(read_only=True))
        self._check_exit_code(diff_command, exit_code=1)
        diff_path = join(f, "diff")
        with open(diff_path) as diff_f:
            diff = diff_f.read()
        for diff_line in DIFF_LINES:
            assert diff_line in diff

    def _make_deterministic(self, node):
        for x in node.findall("testcase"):
            if "time" in x.attrib:
                del x.attrib["time"]

            # Remove contents of stdout/stderr blocks
            for y in x.findall("system-out"):
                y.text = ""
            for y in x.findall("system-err"):
                y.text = ""

        return node
