
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
            self._check_exit_code(["shed_lint"], exit_code=-1)
        with self._isolate_repo("bad_readme_md"):
            self._check_exit_code(["shed_lint"], exit_code=-1)
        with self._isolate_repo("bad_repo_name"):
            self._check_exit_code(["shed_lint"], exit_code=-1)
