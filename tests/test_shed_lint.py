from os.path import (
    join,
    realpath,
)
from shutil import (
    copyfile,
    move,
)

from .test_utils import CliTestCase


class ShedLintTestCase(CliTestCase):
    def test_valid_repos(self):
        with self._isolate_repo("single_tool"):
            self._check_exit_code(["shed_lint", "--skip", "shed_remote_repository_url"])
        with self._isolate_repo("multi_repos_nested"):
            self._check_exit_code(["shed_lint", "--recursive"])
        with self._isolate_repo("package_1"):
            self._check_exit_code(["shed_lint", "--skip", "shed_remote_repository_url"])
        with self._isolate_repo("suite_1"):
            self._check_exit_code(["shed_lint", "--skip", "shed_remote_repository_url"])
        with self._isolate_repo("workflow_1"):
            self._check_exit_code(["shed_lint", "--skip", "shed_remote_repository_url"])

    def test_invalid_repos(self):
        # And now
        with self._isolate_repo("bad_readme_rst"):
            self._check_exit_code(["shed_lint", "--skip", "shed_remote_repository_url"], exit_code=1)
        with self._isolate_repo("bad_readme_md"):
            self._check_exit_code(["shed_lint", "--skip", "shed_remote_repository_url"], exit_code=0)
        with self._isolate_repo("bad_repo_name"):
            self._check_exit_code(["shed_lint", "--skip", "shed_remote_repository_url"], exit_code=1)
        with self._isolate_repo("bad_missing_include"):
            self._check_exit_code(["shed_lint", "--skip", "shed_remote_repository_url"], exit_code=1)
        with self._isolate_repo("bad_missing_tool_deps"):
            self._check_exit_code(["shed_lint", "--skip", "shed_remote_repository_url"], exit_code=1)
        with self._isolate_repo("bad_missing_repo_deps"):
            self._check_exit_code(["shed_lint", "--skip", "shed_remote_repository_url"], exit_code=1)
        with self._isolate_repo("bad_package_category"):
            self._check_exit_code(["shed_lint", "--skip", "shed_remote_repository_url"], exit_code=1)
        with self._isolate_repo("bad_invalid_yaml"):
            self._check_exit_code(["shed_lint", "--skip", "shed_remote_repository_url"], exit_code=254)

    def test_tool_linting(self):
        # Make sure bad_invalid_tool_xml only when used with --tools.
        with self._isolate_repo("bad_invalid_tool_xml"):
            self._check_exit_code(["shed_lint"], exit_code=0)
        with self._isolate_repo("bad_invalid_tool_xml"):
            self._check_exit_code(["shed_lint", "--tools"], exit_code=1)
        with self._isolate_repo("bad_tool_no_citations"):
            self._check_exit_code(["shed_lint", "--tools"], exit_code=1)

    def test_tool_linting_required_files(self):
        # Regression test for https://github.com/galaxyproject/planemo/issues/1646:
        # a sibling file declared in <required_files> must be found in the
        # realized repository even though shed_lint copies files into a temp dir.
        with self._isolate_repo("single_tool_required_files"):
            self._check_exit_code(["shed_lint", "--tools", "--skip", "shed_remote_repository_url"])

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

    def test_fail_fast_on_realization_error(self):
        # A malformed repository definition is rejected immediately.
        with self._isolate() as f:
            for name in ["bad_invalid_yaml", "single_tool_exclude"]:
                self._copy_repo(name, join(f, name))
            r = self._check_exit_code(["shed_lint", "-r", "--fail_fast"], exit_code=1)
            assert isinstance(r.exception, RuntimeError)

    def test_fail_fast_on_lint_failure(self):
        # A normal lint failure should stop both the remaining linters for the
        # repository and traversal of subsequent repositories.
        with self._isolate() as f:
            bad_repo = join(f, "a_bad_repo")
            good_repo = join(f, "b_good_repo")
            self._copy_repo("bad_missing_include", bad_repo)
            self._copy_repo("single_tool", good_repo)

            error_level_result = self._check_exit_code(
                [
                    "shed_lint",
                    "-r",
                    "--fail_fast",
                    "--fail_level",
                    "error",
                    "--skip",
                    "version_bumped",
                ]
            )
            assert ".shed.yml found and appears to be valid YAML." in error_level_result.output
            assert f"Linting repository {realpath(good_repo)}" in error_level_result.output

            r = self._check_exit_code(
                [
                    "shed_lint",
                    "-r",
                    "--fail_fast",
                    "--skip",
                    "version_bumped",
                ],
                exit_code=1,
            )

            assert f"Linting repository {realpath(bad_repo)}" in r.output
            assert "Failed to expand inclusions" in r.output
            assert ".shed.yml found and appears to be valid YAML." not in r.output
            assert f"Linting repository {realpath(good_repo)}" not in r.output
            assert not isinstance(r.exception, RuntimeError)
            assert "Traceback" not in r.output

    def test_fail_fast_while_linting_tools(self):
        # A repository can contain multiple tools, so stopping repository
        # traversal alone is not sufficient to implement --fail_fast.
        with self._isolate_repo("single_tool") as f:
            first_bad_tool = join(f, "first_bad_tool.xml")
            second_bad_tool = join(f, "second_bad_tool.xml")
            move(join(f, "cat.xml"), first_bad_tool)
            copyfile(first_bad_tool, second_bad_tool)

            r = self._check_exit_code(
                ["shed_lint", "--tools", "--fail_fast", "--skip", "version_bumped"],
                exit_code=1,
            )

            assert r.output.count("+Linting tool ") == 1
            assert "No citations found" in r.output
            assert "Traceback" not in r.output

    def test_ensure_metadata(self):
        with self._isolate_repo("single_tool"):
            self._check_exit_code(["shed_lint", "--skip", "shed_remote_repository_url"])
        with self._isolate_repo("single_tool_exclude"):
            self._check_exit_code(
                ["shed_lint", "--skip", "shed_remote_repository_url", "--ensure_metadata"], exit_code=1
            )
