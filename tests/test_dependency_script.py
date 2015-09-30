import os
import glob

from .test_utils import (
    CliTestCase,
    TEST_REPOS_DIR,
)


class LintTestCase(CliTestCase):

    def test_empty_file(self):
        # Using empty file repos/fastqc/tool_dependencies.xml
        ds_cmd = ["dependency_script",
                  os.path.join(TEST_REPOS_DIR, "fastqc")]
        self._check_exit_code(ds_cmd, 1)

    def test_samtools_example(self):
        ds_cmd = ["dependency_script",
                  os.path.join(TEST_REPOS_DIR, "package_1")]
        self._check_exit_code(ds_cmd, 0)

    def test_cd_hit_auxtools(self):
        ds_cmd = ["dependency_script",
                  os.path.join(TEST_REPOS_DIR, "bad_package_name")]
        self._check_exit_code(ds_cmd, 0)

    def test_good_examples(self):
        ds_cmd = ["dependency_script",
                  os.path.join(TEST_REPOS_DIR, "package_1"),
                  os.path.join(TEST_REPOS_DIR, "bad_package_category/"),
                  os.path.join(TEST_REPOS_DIR, "bad_package_name")]
        self._check_exit_code(ds_cmd, 0)

    def test_repos_recurse(self):
        # At least one will fail
        ds_cmd = ["dependency_script", "-r", TEST_REPOS_DIR]
        self._check_exit_code(ds_cmd, 1)

    def test_repos_recurse_fast(self):
        # At least one will fail
        ds_cmd = ["dependency_script", "--fail_fast", "-r", TEST_REPOS_DIR]
        self._check_exit_code(ds_cmd, 1)

    def test_repos_list(self):
        # At least one will fail
        ds_cmd = ["dependency_script"]  + list(glob.glob("%s/*" % TEST_REPOS_DIR))
        self._check_exit_code(ds_cmd, 1)

    def test_repos_list_fast(self):
        # At least one will fail
        ds_cmd = ["dependency_script", "--fail_fast"]  + list(glob.glob("%s/*" % TEST_REPOS_DIR))
        self._check_exit_code(ds_cmd, 1)
