from os.path import join

from .test_utils import CliTestCase


class ShedLintTestCase(CliTestCase):
    def test_valid_repos(self):
        with self._isolate_repo("single_tool"):
            self._check_exit_code(["shed_lint"])
        with self._isolate_repo("multi_repos_nested"):
            self._check_exit_code(["shed_lint", "--recursive"])
        with self._isolate_repo("package_1"):
            self._check_exit_code(["shed_lint"])
        with self._isolate_repo("suite_1"):
            self._check_exit_code(["shed_lint"])
        with self._isolate_repo("workflow_1"):
            self._check_exit_code(["shed_lint"])

    def test_invalid_repos(self):
        # And now
        with self._isolate_repo("bad_readme_rst"):
            self._check_exit_code(["shed_lint"], exit_code=1)
        with self._isolate_repo("bad_readme_md"):
            self._check_exit_code(["shed_lint"], exit_code=1)
        with self._isolate_repo("bad_repo_name"):
            self._check_exit_code(["shed_lint"], exit_code=1)
        with self._isolate_repo("bad_missing_include"):
            self._check_exit_code(["shed_lint"], exit_code=1)
        with self._isolate_repo("bad_missing_tool_deps"):
            self._check_exit_code(["shed_lint"], exit_code=1)
        with self._isolate_repo("bad_missing_repo_deps"):
            self._check_exit_code(["shed_lint"], exit_code=1)
        with self._isolate_repo("bad_package_category"):
            self._check_exit_code(["shed_lint"], exit_code=1)
        with self._isolate_repo("bad_invalid_yaml"):
            self._check_exit_code(["shed_lint"], exit_code=254)

    def test_tool_linting(self):
        # Make sure bad_invalid_tool_xml only when used with --tools.
        with self._isolate_repo("bad_invalid_tool_xml"):
            self._check_exit_code(["shed_lint"], exit_code=0)
        with self._isolate_repo("bad_invalid_tool_xml"):
            self._check_exit_code(["shed_lint", "--tools"], exit_code=1)
        with self._isolate_repo("bad_tool_no_citations"):
            self._check_exit_code(["shed_lint", "--tools"], exit_code=1)

    def test_invalid_nested(self):
        # Created a nested repository with one good and one
        # invalid repository and make sure it runs and produces
        # a 254 (it ran to completion but one or more things failed
        # )
        with self._isolate() as f:
            for name in ["bad_invalid_yaml", "single_tool_exclude"]:
                self._copy_repo(name, join(f, name))
                self._copy_repo(name, join(f, name))
            self._check_exit_code(["shed_lint", "-r"], exit_code=254)

    def test_fail_fast(self):
        # Created a nested repository with one good and one
        # invalid repository and make sure it exits immediately with 1.
        with self._isolate() as f:
            for name in ["bad_invalid_yaml", "single_tool_exclude"]:
                self._copy_repo(name, join(f, name))
                self._copy_repo(name, join(f, name))
            r = self._check_exit_code(["shed_lint", "-r", "--fail_fast"], exit_code=1)
            assert isinstance(r.exception, RuntimeError)

    def test_ensure_metadata(self):
        with self._isolate_repo("single_tool"):
            self._check_exit_code(["shed_lint"])
        with self._isolate_repo("single_tool_exclude"):
            self._check_exit_code(["shed_lint", "--ensure_metadata"], exit_code=1)
