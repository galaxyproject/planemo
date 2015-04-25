""" Integration tests for shed_upload, shed_download, and shed_create
commands.
"""
from os.path import exists, join
import os
import tarfile
import shutil

from .test_utils import CliShedTestCase


class ShedUploadTestCase(CliShedTestCase):

    def test_tar_single(self):
        with self._isolate_repo("single_tool") as f:
            upload_command = ["shed_upload", "--tar_only"]
            upload_command.extend(self._shed_args())
            self._check_exit_code(upload_command)
            assert exists(join(f, "shed_upload.tar.gz"))

    def test_upload_not_exists(self):
        with self._isolate_repo("single_tool"):
            upload_command = ["shed_upload"]
            upload_command.extend(self._shed_args())
            self._check_exit_code(upload_command, exit_code=-1)

    def test_upload_with_force_create(self):
        with self._isolate_repo("single_tool") as f:
            upload_command = ["shed_upload", "--force_repository_creation"]
            upload_command.extend(self._shed_args())
            self._check_exit_code(upload_command)
            self._verify_single_uploaded(f)

    def test_create_and_upload(self):
        with self._isolate_repo("single_tool") as f:
            create_command = ["shed_create"]
            create_command.extend(self._shed_args())
            self._check_exit_code(create_command)
            upload_command = ["shed_upload"]
            upload_command.extend(self._shed_args())
            self._check_exit_code(upload_command)
            self._verify_single_uploaded(f)

    def test_cannont_recreate(self):
        with self._isolate_repo("single_tool"):
            create_command = ["shed_create"]
            create_command.extend(self._shed_args())
            self._check_exit_code(create_command)
            self._check_exit_code(create_command, exit_code=-1)

    def test_upload_recusrive(self):
        with self._isolate_repo("multi_repos_nested") as f:
            upload_command = [
                "shed_upload", "-r", "--force_repository_creation"
            ]
            upload_command.extend(self._shed_args())
            self._check_exit_code(upload_command)
            self._verify_upload(f, ["cat1.xml", "macros.xml"], ["cat1"])
            self._verify_upload(f, ["cat2.xml", "macros.xml"], ["cat2"])

    def test_upload_filters_invalid_suite(self):
        with self._isolate_repo("suite_1") as f:
            # No .shed.yml, make sure to test it can infer type
            # from passed in --name.
            upload_command = ["shed_upload", "--tar_only",
                              "--name", "suite_1"]
            upload_command.extend(self._shed_args())
            self._check_exit_code(upload_command)
            target = self._untar(f, "shed_upload.tar.gz")
            # Only one file was in archive
            assert exists(join(target, "repository_dependencies.xml"))
            # this got filtered
            assert not exists(join(target, "README.rst"))

    def test_upload_filters_ignore(self):
        with self._isolate_repo("single_tool_exclude") as f:
            upload_command = ["shed_upload", "--force_repository_creation"]
            upload_command.extend(self._shed_args())
            self._check_exit_code(upload_command)
            target = self._verify_upload(f)
            assert not exists(join(target, "related_file"))

    def test_tar_with_symlinks(self):
        with self._isolate_repo("multi_repos_nested") as f:
            upload_command = ["shed_upload", "--force_repository_creation"]
            upload_command.append("cat2")
            upload_command.extend(self._shed_args())
            self._check_exit_code(upload_command)
            target = self._verify_upload(f, ["macros.xml"], ["cat2"])
            with open(join(target, "macros.xml"), "r") as macro_f:
                macro_contents = macro_f.read()
                assert macro_contents.startswith("<macros>")

    def test_upload_filters_git(self):
        with self._isolate_repo("single_tool") as f:
            mock_git_dir = join(f, ".git")
            os.makedirs(mock_git_dir)
            index_path = os.path.join(mock_git_dir, "index_file")
            with open(index_path, "w") as index_f:
                index_f.write("test")
            upload_command = ["shed_upload", "--force_repository_creation"]
            upload_command.extend(self._shed_args())
            self._check_exit_code(upload_command)
            target = self._verify_upload(f)
            assert not exists(join(target, ".git"))

    def test_upload_filters_invalid_package(self):
        with self._isolate_repo("package_1") as f:
            upload_command = ["shed_upload", "--tar_only"]
            upload_command.extend(self._shed_args())
            self._check_exit_code(upload_command)
            target = self._untar(f, "shed_upload.tar.gz")
            # Only one file was in archive
            assert exists(join(target, "tool_dependencies.xml"))
            # this got filtered
            assert not exists(join(target, "README.rst"))
            # .shed.yml always gets filtered
            assert not exists(join(target, ".shed.yml"))

    def test_upload_not_filters_unrestricted(self):
        with self._isolate_repo("workflow_1") as f:
            upload_command = ["shed_upload", "--tar_only"]
            upload_command.extend(self._shed_args())
            self._check_exit_code(upload_command)
            target = self._untar(f, "shed_upload.tar.gz")
            # Only one file was in archive
            assert exists(join(target, "repository_dependencies.xml"))
            assert exists(join(target, "README.rst"))

    def test_upload_expansion_configured(self):
        with self._isolate_repo("multi_repos_flat_configured") as f:
            self._verify_expansion(f)

    def test_upload_expansion_flagged(self):
        with self._isolate_repo("multi_repos_flat_flag") as f:
            self._verify_expansion(f)

    def _verify_expansion(self, f):
        upload_command = ["shed_upload", "--tar_only"]
        upload_command.extend(self._shed_args())
        self._check_exit_code(upload_command)
        self._check_tar(
            f, "shed_upload_cs_cat1.tar.gz",
            contains=[
                "cat1.xml",
                "macros.xml",
                "test-data/1.bed"
            ],
            not_contains=["cat2.xml"]
        )
        self._check_tar(
            f, "shed_upload_cs_cat2.tar.gz",
            contains=[
                "cat2.xml",
                "macros.xml",
                "test-data/1.bed"
            ],
            not_contains=["cat1.xml"]
        )

    def _verify_single_uploaded(self, f):
        self._verify_upload(f, ["cat.xml", "related_file", "test-data/1.bed"])

    def _verify_upload(self, f, download_files=[], download_args=[]):
        target = self._download_repo(f, download_args)
        for download_file in download_files:
            assert exists(join(target, download_file)), download_file
        return target

    def _check_tar(self, f, tar_path, contains=[], not_contains=[]):
        tar_path = join(f, tar_path)
        assert exists(tar_path)
        target = self._untar(f, tar_path)
        for path in contains:
            assert exists(join(target, path))
        for path in not_contains:
            assert not exists(join(target, path))

    def _download_repo(self, f, download_args=[]):
        download_command = ["shed_download"]
        download_command.extend(download_args)
        download_command.extend(self._shed_args(read_only=True))
        self._check_exit_code(download_command)
        download = join(f, "shed_download.tar.gz")
        assert exists(download)
        return self._untar(f, "shed_download.tar.gz")

    def _untar(self, f, path):
        target = join(f, "download")
        if exists(target):
            shutil.rmtree(target)
        os.makedirs(target)
        tar = tarfile.open(path, "r:gz")
        tar.extractall(path=target)
        tar.close()
        return target
