
from .test_utils import CliTestCase


class ShedLineTestCase(CliTestCase):

    def test_valid_repo(self):
        with self._isolate_repo("single_tool"):
            self._check_exit_code(["shed_lint"])
