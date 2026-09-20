from os.path import join

from .test_utils import (
    assert_exists,
    CliShedTestCase,
)


class ShedUploadTestCase(CliShedTestCase):
    def test_download_expansion_tars(self):
        with self._isolate_repo("multi_repos_flat_flag") as f:
            upload_command = ["shed_upload", "--force_repository_creation"]
            upload_command.extend(self._shed_args())
            self._check_exit_code(upload_command)
            download_command = [
                "shed_download",
            ]
            download_command.extend(self._shed_args(read_only=True))
            self._check_exit_code(download_command)
            assert_exists(join(f, "shed_download_cs_cat1.tar.gz"))
            assert_exists(join(f, "shed_download_cs_cat2.tar.gz"))

    def test_download_expansion_dir(self):
        with self._isolate_repo("multi_repos_flat_flag") as f:
            upload_command = ["shed_upload", "--force_repository_creation"]
            upload_command.extend(self._shed_args())
            self._check_exit_code(upload_command)
            download_command = ["shed_download", "--destination", "shed_d"]
            download_command.extend(self._shed_args(read_only=True))
            self._check_exit_code(download_command)
            assert_exists(join(f, "shed_d_cs_cat1", "cat1.xml"))
            assert_exists(join(f, "shed_d_cs_cat2", "cat2.xml"))
