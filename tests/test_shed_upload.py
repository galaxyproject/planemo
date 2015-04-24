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
            self._check_exit_code(create_command, exit_code=1)

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

    def _verify_single_uploaded(self, f):
        self._verify_upload(f, ["cat.xml"])

    def _verify_upload(self, f, download_files=[], download_args=[]):
        target = self._download_repo(f, download_args)
        for download_file in download_files:
            assert exists(join(target, download_file))

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
