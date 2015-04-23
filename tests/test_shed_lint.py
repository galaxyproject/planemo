
from .test_utils import CliTestCase


class ShedLineTestCase(CliTestCase):

    def test_valid_repo(self):
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
