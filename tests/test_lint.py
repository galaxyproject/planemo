import glob
import os

from .test_utils import (
    CliTestCase,
    TEST_REPOS_DIR,
    TEST_TOOLS_DIR,
)


class LintTestCase(CliTestCase):

    def test_ok_tools(self):
        ok_tools = glob.glob("%s/ok_*" % TEST_TOOLS_DIR)
        for ok_tool in ok_tools:
            lint_cmd = ["lint", "--urls", ok_tool]
            self._check_exit_code(lint_cmd)

    def test_fail_tools(self):
        fail_tools = glob.glob("%s/fail_*" % TEST_TOOLS_DIR)
        for fail_tool in fail_tools:
            lint_cmd = ["lint", fail_tool]
            self._check_exit_code(lint_cmd, exit_code=1)

    def test_lint_default(self):
        with self._isolate_repo("single_tool"):
            self._check_exit_code(["lint", "--skip", "citations"])
        with self._isolate_repo("single_tool"):
            self._check_exit_code(["lint"], exit_code=1)

    def test_lint_multiple(self):
        names = ["fail_citation.xml", "fail_order.xml"]
        paths = list(map(lambda p: os.path.join(TEST_TOOLS_DIR, p), names))
        self._check_exit_code(["lint"] + paths, exit_code=1)
        self._check_exit_code(
            ["lint", "--skip", "citations,xml_order"] + paths,
            exit_code=0
        )

    def test_skips(self):
        fail_citation = os.path.join(TEST_TOOLS_DIR, "fail_citation.xml")
        lint_cmd = ["lint", fail_citation]
        self._check_exit_code(lint_cmd, exit_code=1)

        lint_cmd = ["lint", "--skip", "citations", fail_citation]
        self._check_exit_code(lint_cmd, exit_code=0)

        # Check string splitting and stuff.
        lint_cmd = ["lint", "--skip", "xml_order, citations", fail_citation]
        self._check_exit_code(lint_cmd, exit_code=0)

    def test_recursive(self):
        nested_dir = os.path.join(TEST_REPOS_DIR, "multi_repos_nested")

        # Fails to find any tools without -r.
        lint_cmd = ["lint", "--skip", "citations", nested_dir]
        self._check_exit_code(lint_cmd, exit_code=2)

        # Works with -r.
        lint_cmd = ["lint", "--skip", "citations", "-r", nested_dir]
        self._check_exit_code(lint_cmd, exit_code=0)

    def test_empty_cdata(self):
        empty_cdata = os.path.join(TEST_TOOLS_DIR, "empty_cdata.xml")
        lint_cmd = ["lint", "--skip", "citations,help", empty_cdata]
        self._check_exit_code(lint_cmd, exit_code=0)

    def test_lint_doi(self):
        fail_doi = os.path.join(TEST_TOOLS_DIR, "invalid_doi.xml")
        lint_cmd = ["lint", "--doi", fail_doi]
        self._check_exit_code(lint_cmd, exit_code=1)
