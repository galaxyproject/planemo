import os
import shutil
import tarfile

from .test_utils import (
    TempDirectoryTestCase,
    TEST_REPOS_DIR,
)

from planemo import shed


class ShedExpansionTestCase(TempDirectoryTestCase):

    def setUp(self):
        super(ShedExpansionTestCase, self).setUp()
        self.target = os.path.join(self.temp_directory, "target")
        os.makedirs(self.target)

    def test_tar_simple(self):
        single_tool = os.path.join(TEST_REPOS_DIR, "single_tool")
        self._tar_and_expand(single_tool)
        cat_path = self._assert_target_exists("cat.xml")
        assert open(cat_path, "r").read().startswith("<tool")
        self._assert_target_exists("related_file")
        self._assert_target_exists("test-data", "1.bed")

    def test_tar_ignore(self):
        single_tool = os.path.join(TEST_REPOS_DIR, "single_tool_exclude")
        self._tar_and_expand(single_tool)
        self._assert_target_not_exists("related_file")

    def test_tar_with_symlinks(self):
        cat2_dir = os.path.join(TEST_REPOS_DIR, "multi_repos_nested", "cat2")
        self._tar_and_expand(cat2_dir)
        macros_path = self._assert_target_exists("macros.xml")
        assert open(macros_path, "r").read().startswith("<macros>")

    def test_git_ignored(self):
        source = os.path.join(self.temp_directory, "src")
        single_tool = os.path.join(TEST_REPOS_DIR, "single_tool")
        shutil.copytree(single_tool, source)
        mock_git_dir = os.path.join(source, ".git")
        os.makedirs(mock_git_dir)
        open(os.path.join(mock_git_dir, "index_file"), "w").write("test")
        self._tar_and_expand(source)
        self._assert_target_not_exists(".git")

    def _assert_target_exists(self, *args):
        path = self._target_path(*args)
        assert os.path.exists(path)
        return path

    def _assert_target_not_exists(self, *args):
        assert not os.path.exists(self._target_path(*args))

    def _target_path(self, *args):
        return os.path.join(self.target, *args)

    def _tar_and_expand(self, test_directory):
        tar_path = shed.build_tarball(test_directory)
        tar_file = tarfile.open(name=tar_path, mode="r")
        tar_file.extractall(path=self.target)
