import glob
import os
import shutil
import tempfile

from .test_utils import (
    CliTestCase,
    skip_if_environ,
    TEST_REPOS_DIR,
)


class DependencyScriptTestCase(CliTestCase):

    @classmethod
    def setUpClass(clazz):
        # This is a bit of a hack, but having failed to get the test
        # script to see $INSTALL_DIR and $DOWNLOAD_CACHE set via
        # .travis.yml this is an alternative.
        temp_dir = tempfile.mkdtemp()
        clazz.install_dir = os.path.join(temp_dir, "install_dir")
        clazz.download_cache = os.path.join(temp_dir, "download_cache")

        for env_var, default in [
                ("INSTALL_DIR", clazz.install_dir),
                ("DOWNLOAD_CACHE", clazz.download_cache),
        ]:
            if env_var not in os.environ:
                os.environ[env_var] = default
                os.makedirs(default)
            value = os.environ[env_var]
            assert os.path.isdir(value), value

    @classmethod
    def tearDownClass(clazz):
        shutil.rmtree(clazz.install_dir)
        shutil.rmtree(clazz.download_cache)

    def setUp(self):
        # Pass to base-class
        CliTestCase.setUp(self)
        # TODO: Switch this to setUpClass once drop Python 2.6
        repo_list = []
        for x in glob.glob("%s/*" % TEST_REPOS_DIR):
            if os.path.isdir(x):
                repo_list.append(x)
        self.repo_list = repo_list

    def test_empty_file(self):
        with self._isolate():
            # Using empty file repos/fastqc/tool_dependencies.xml
            ds_cmd = ["dependency_script",
                      os.path.join(TEST_REPOS_DIR, "fastqc")]
            self._check_exit_code(ds_cmd, 1)

    @skip_if_environ("PLANEMO_SKIP_SLOW_TESTS")
    def test_samtools_example(self):
        with self._isolate():
            # Also checking the --download_cache option
            ds_cmd = ["dependency_script",
                      "--download_cache", os.environ["DOWNLOAD_CACHE"],
                      os.path.join(TEST_REPOS_DIR, "package_1")]
            self._check_exit_code(ds_cmd, 0)

    @skip_if_environ("PLANEMO_SKIP_SLOW_TESTS")
    def test_cd_hit_auxtools(self):
        with self._isolate():
            ds_cmd = ["dependency_script",
                      os.path.join(TEST_REPOS_DIR, "bad_repo_name")]
            self._check_exit_code(ds_cmd, 0)

    @skip_if_environ("PLANEMO_SKIP_SLOW_TESTS")
    def test_good_examples(self):
        with self._isolate():
            ds_cmd = ["dependency_script",
                      os.path.join(TEST_REPOS_DIR, "package_1"),
                      os.path.join(TEST_REPOS_DIR, "bad_package_category/"),
                      os.path.join(TEST_REPOS_DIR, "bad_repo_name")]
            self._check_exit_code(ds_cmd, 0)

    @skip_if_environ("PLANEMO_SKIP_SLOW_TESTS")
    def test_repos_recurse(self):
        with self._isolate():
            # At least one will fail
            ds_cmd = ["dependency_script", "-r", TEST_REPOS_DIR]
            self._check_exit_code(ds_cmd, 1)

    @skip_if_environ("PLANEMO_SKIP_SLOW_TESTS")
    def test_repos_recurse_fast(self):
        with self._isolate():
            # At least one will fail
            ds_cmd = ["dependency_script", "--fail_fast", "-r", TEST_REPOS_DIR]
            self._check_exit_code(ds_cmd, 1)

    @skip_if_environ("PLANEMO_SKIP_SLOW_TESTS")
    def test_repos_list(self):
        with self._isolate():
            # At least one will fail
            ds_cmd = ["dependency_script"] + self.repo_list
            self._check_exit_code(ds_cmd, 1)

    @skip_if_environ("PLANEMO_SKIP_SLOW_TESTS")
    def test_repos_list_fast(self):
        with self._isolate():
            # At least one will fail
            ds_cmd = ["dependency_script", "--fail_fast"] + self.repo_list
            self._check_exit_code(ds_cmd, 1)
