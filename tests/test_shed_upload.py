""" Integration tests for shed_upload, shed_download, and shed_create
commands.
"""
import os
import tarfile

from .test_utils import CliShedTestCase


class ShedUploadTestCase(CliShedTestCase):

    def test_tar_single(self):
        with self._isolate_repo("single_tool") as f:
            upload_command = ["shed_upload", "--tar_only"]
            upload_command.extend(self._shed_args())
            self._check_exit_code(upload_command)
            assert os.path.exists(os.path.join(f, "shed_upload.tar.gz"))

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

    def _verify_single_uploaded(self, f):
        download_command = ["shed_download"]
        download_command.extend(self._shed_args(read_only=True))
        self._check_exit_code(download_command)
        download = os.path.join(f, "shed_download.tar.gz")
        assert os.path.exists(download)
        target = os.path.join(f, "download")
        os.makedirs(target)
        tar = tarfile.open(download, "r:gz")
        tar.extractall(path=target)
        tar.close()
        assert os.path.exists(os.path.join(target, "cat.xml"))
