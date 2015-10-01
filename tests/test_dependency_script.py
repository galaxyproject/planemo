import os
import glob
import tempfile

from .test_utils import (
    CliTestCase,
    TEST_REPOS_DIR,
)


# This is a bit of a hack, but having failed to get the test
# script to see $INSTALL_DIR and $DOWNLOAD_CACHE set via
# .travis.yml this is an alternative.
temp_dir = tempfile.mkdtemp()
for env_var, default in [
        ("INSTALL_DIR", os.path.join(temp_dir, "install_dir")),
        ("DOWNLOAD_CACHE", os.path.join(temp_dir, "download_cache")),
        ]:
    if env_var not in os.environ:
        os.environ[env_var] = default
        os.makedirs(default)
    value = os.environ[env_var]
    print("Using $%s=%s" % (env_var, value))
    assert os.path.isdir(value), value
# TODO - Remove temp_dir


class DependencyScriptTestCase(CliTestCase):

    def test_empty_file(self):
        # Using empty file repos/fastqc/tool_dependencies.xml
        ds_cmd = ["dependency_script",
                  os.path.join(TEST_REPOS_DIR, "fastqc")]
        self._check_exit_code(ds_cmd, 1)

    def test_samtools_example(self):
        # Also checking the --download_cache option
        ds_cmd = ["dependency_script",
                  "--download_cache", os.environ["DOWNLOAD_CACHE"],
                  os.path.join(TEST_REPOS_DIR, "package_1")]
        self._check_exit_code(ds_cmd, 0)

    def test_cd_hit_auxtools(self):
        ds_cmd = ["dependency_script",
                  os.path.join(TEST_REPOS_DIR, "bad_repo_name")]
        self._check_exit_code(ds_cmd, 0)

    def test_good_examples(self):
        ds_cmd = ["dependency_script",
                  os.path.join(TEST_REPOS_DIR, "package_1"),
                  os.path.join(TEST_REPOS_DIR, "bad_package_category/"),
                  os.path.join(TEST_REPOS_DIR, "bad_repo_name")]
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
        ds_cmd = ["dependency_script"] + list(glob.glob("%s/*" % TEST_REPOS_DIR))
        self._check_exit_code(ds_cmd, 1)

    def test_repos_list_fast(self):
        # At least one will fail
        ds_cmd = ["dependency_script", "--fail_fast"] + list(glob.glob("%s/*" % TEST_REPOS_DIR))
        self._check_exit_code(ds_cmd, 1)
