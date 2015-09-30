import os
import glob

from .test_utils import (
    CliTestCase,
    TEST_REPOS_DIR,
)


class LintTestCase(CliTestCase):

    def test_repos_recurse(self):
        ds_cmd = ["dependency_script", "-r", TEST_REPOS_DIR]
        self._check_exit_code(ds_cmd)

    def test_repos_individually(self):
        for repo in glob.glob("%s/*" % TEST_REPOS_DIR):
            ds_cmd = ["dependency_script", repo]
            self._check_exit_code(ds_cmd)
