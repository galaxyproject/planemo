"""Integration tests for shed contents commands.

Specifically, tests for shed_upload, shed_download, and shed_create.
commands.
"""
import contextlib
import os
import shutil
import tarfile
from os.path import (
    exists,
    join,
)

from galaxy.util import unicodify

from planemo.io import shell
from .test_utils import (
    assert_exists,
    CliShedTestCase,
    modify_environ,
    TEST_REPOS_DIR,
)


class ShedUploadTestCase(CliShedTestCase):
    def test_tar_single(self):
        with self._isolate_repo("single_tool") as f:
            upload_command = ["shed_upload", "--tar_only"]
            upload_command.extend(self._shed_args())
            self._check_exit_code(upload_command)
            assert_exists(join(f, "shed_upload.tar.gz"))

    def test_upload_not_exists(self):
        with self._isolate_repo("single_tool"):
            upload_command = ["shed_upload"]
            upload_command.extend(self._shed_args())
            self._check_exit_code(upload_command, exit_code=2)

    def test_update_not_exists(self):
        with self._isolate_repo("single_tool"):
            upload_command = ["shed_update"]
            upload_command.extend(self._shed_args())
            self._check_exit_code(upload_command, exit_code=2)

    def test_update_not_exists_update_only(self):
        with self._isolate_repo("single_tool"):
            upload_command = ["shed_update", "--skip_upload"]
            upload_command.extend(self._shed_args())
            self._check_exit_code(upload_command, exit_code=2)

    def test_update_with_check_diff(self):
        with self._isolate_repo("single_tool") as f:
            self._shed_create()

            self._assert_shed_diff(diff=0)

            upload_command = ["shed_update", "--force_repository_creation", "--check_diff"]
            upload_command.extend(self._shed_args())
            self._check_exit_code(upload_command)

            # First time no difference.
            r = self._check_exit_code(upload_command)
            assert "not different, skipping upload." in r.output

            # Modify a file so there is a difference.
            with open(join(f, "related_file"), "w") as rf:
                rf.write("new_contents")

            self._assert_shed_diff(diff=1)

            # No assert there is no difference again.
            r = self._check_exit_code(upload_command)
            assert "not different, skipping upload." not in r.output

            self._assert_shed_diff(diff=0)

    def test_update_with_check_diff_package(self):
        with self._isolate_repo("package_1") as f:
            self._shed_create()

            self._assert_shed_diff(diff=0)
            upload_command = ["shed_update", "--force_repository_creation", "--check_diff"]
            upload_command.extend(self._shed_args())
            self._check_exit_code(upload_command)

            # First time no difference.
            r = self._check_exit_code(upload_command)
            assert "not different, skipping upload." in r.output

            update_package_1(f)
            self._assert_shed_diff(diff=1)

            # No assert there is no difference again.
            r = self._check_exit_code(upload_command)
            assert "not different, skipping upload." not in r.output

            self._assert_shed_diff(diff=0)

    def test_update_with_force_create_metadata_only(self):
        with self._isolate_repo("single_tool") as f:
            upload_command = ["shed_update", "--force_repository_creation", "--skip_upload"]
            upload_command.extend(self._shed_args())
            self._check_exit_code(upload_command)
            self._verify_empty_repository(f)

    def test_update_with_force_create(self):
        with self._isolate_repo("single_tool") as f:
            upload_command = ["shed_update", "--force_repository_creation"]
            upload_command.extend(self._shed_args())
            self._check_exit_code(upload_command)
            self._verify_single_uploaded(f)

    def test_tar_from_git(self):
        with self._isolate() as f:
            with self._git_configured():
                dest = join(f, "single_tool")
                self._copy_repo("single_tool", dest)
                shell(" && ".join(["cd %s" % dest, "git init", "git add .", "git commit -m 'initial commit'"]))
                upload_command = ["shed_update", "--force_repository_creation", "git+single_tool/.git"]
                upload_command.extend(self._shed_args())
                self._check_exit_code(upload_command)
                self._verify_single_uploaded(f, ["single_tool"])

    def test_upload_from_git(self):
        with self._isolate() as f:
            with self._git_configured():
                dest = join(f, "single_tool")
                self._copy_repo("single_tool", dest)
                shell(" && ".join(["cd %s" % dest, "git init", "git add .", "git commit -m 'initial commit'"]))
                upload_command = ["shed_update", "--force_repository_creation", "git+single_tool/.git"]
                upload_command.extend(self._shed_args())
                self._check_exit_code(upload_command)
                self._verify_single_uploaded(f, ["single_tool"])

    @contextlib.contextmanager
    def _git_configured(self):
        with modify_environ(
            {
                "GIT_AUTHOR_NAME": "planemo developer",
                "GIT_COMMITTER_NAME": "planemo developer",
                "EMAIL": "planemo@galaxyproject.org",
                "GIT_AUTHOR_EMAIL": "planemo@galaxyproject.org",
                "GIT_COMMITTER_EMAIL": "planemo@galaxyproject.org",
            }
        ):
            yield

    def test_create_and_upload(self):
        with self._isolate_repo("single_tool") as f:
            create_command = ["shed_create"]
            create_command.extend(self._shed_args())
            self._check_exit_code(create_command)
            upload_command = ["shed_upload"]
            upload_command.extend(self._shed_args())
            self._check_exit_code(upload_command)
            self._verify_single_uploaded(f)

    def test_create_with_upload(self):
        with self._isolate_repo("single_tool") as f:
            create_command = ["shed_create"]
            create_command.extend(self._shed_args())
            self._check_exit_code(create_command)
            self._verify_single_uploaded(f)

    def test_cannont_recreate(self):
        with self._isolate_repo("single_tool"):
            create_command = ["shed_create"]
            create_command.extend(self._shed_args())
            self._check_exit_code(create_command)
            self._check_exit_code(create_command, exit_code=1)

    def test_cannot_upload_missing_include(self):
        with self._isolate_repo("bad_missing_include"):
            upload_command = ["shed_upload", "--tar_only"]
            upload_command.extend(self._shed_args())
            self._check_exit_code(upload_command, exit_code=254)

    def test_upload_recusrive(self):
        with self._isolate_repo("multi_repos_nested") as f:
            upload_command = ["shed_update", "-r", "--force_repository_creation"]
            upload_command.extend(self._shed_args())
            self._check_exit_code(upload_command)
            self._verify_upload(f, ["cat1.xml", "macros.xml"], ["cat1"])
            self._verify_upload(f, ["cat2.xml", "macros.xml"], ["cat2"])

    def test_upload_filters_invalid_suite(self):
        with self._isolate_repo("suite_1") as f:
            # No .shed.yml, make sure to test it can infer type
            # from passed in --name.
            upload_command = ["shed_upload", "--tar_only", "--name", "suite_1"]
            upload_command.extend(self._shed_args())
            self._check_exit_code(upload_command)
            target = self._untar(f, "shed_upload.tar.gz")
            # Only one file was in archive
            assert_exists(join(target, "repository_dependencies.xml"))
            # this got filtered
            assert not exists(join(target, "README.rst"))

    def test_upload_suite_auto(self):
        with self._isolate_repo("suite_auto") as f:
            upload_command = ["shed_upload", "--tar_only"]
            upload_command.extend(self._shed_args())
            self._check_exit_code(upload_command)
            target = self._untar(f, "shed_upload.tar.gz")
            # Only one file was in archive
            assert_exists(join(target, "repository_dependencies.xml"))

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
            with open(join(target, "macros.xml")) as macro_f:
                macro_contents = macro_f.read()
                assert macro_contents.startswith("<macros>")

    def test_upload_filters_git(self):
        with self._isolate_repo("single_tool") as f:
            mock_git_dir = join(f, ".git")
            os.makedirs(mock_git_dir)
            index_path = join(mock_git_dir, "index_file")
            with open(index_path, "w") as index_f:
                index_f.write("test")
            with open(join(f, "related_file~"), "w") as tilde_f:
                tilde_f.write("backup!")
            upload_command = ["shed_upload", "--force_repository_creation"]
            upload_command.extend(self._shed_args())
            self._check_exit_code(upload_command)
            target = self._verify_upload(f)
            assert not exists(join(target, ".git"))
            assert not exists(join(target, "related_file~"))

    def test_upload_filters_invalid_package(self):
        with self._isolate_repo("package_1") as f:
            upload_command = ["shed_upload", "--tar_only"]
            upload_command.extend(self._shed_args())
            self._check_exit_code(upload_command)
            target = self._untar(f, "shed_upload.tar.gz")
            # Only one file was in archive
            assert_exists(join(target, "tool_dependencies.xml"))
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
            assert_exists(join(target, "repository_dependencies.xml"))
            assert_exists(join(target, "README.rst"))

    def test_upload_expansion_configured(self):
        with self._isolate_repo("multi_repos_flat_configured") as f:
            self._verify_expansion(f)

    def test_upload_expansion_flagged(self):
        with self._isolate_repo("multi_repos_flat_flag") as f:
            self._verify_expansion(f)

    def test_upload_expansion_configured_extras(self):
        with self._isolate() as f:
            repo = join(f, "repo")
            self._copy_repo("multi_repos_flat_configured_complex", repo)
            self._copy_repo("shared_files", join(f, "shared_files"))
            self._verify_expansion(f, "repo")
            for tool_id in ["cat1", "cat2"]:
                self._check_tar(
                    f,
                    "shed_upload_cs_%s.tar.gz" % tool_id,
                    contains=[
                        "CITATION.txt",
                        "test-data/extra_test_file.txt",
                    ],
                )

    def test_upload_expansion_suite(self):
        with self._isolate_repo("multi_repos_flat_flag_suite") as f:
            self._verify_expansion(f)
            target = self._check_tar(
                f,
                "shed_upload_suite_cat.tar.gz",
                contains=[
                    "repository_dependencies.xml",
                ],
                not_contains=["macros.xml"],
            )
            with open(join(target, "repository_dependencies.xml")) as f:
                repo_xml = f.read()
                assert 'owner="devteam" name="cat_legacy"' in repo_xml
                assert 'owner="iuc" name="cs-cat2"' in repo_xml

    def test_upload_with_double_dot(self):
        with self._isolate() as f:
            self._copy_repo("up_root/", join(f, "up_root/"))
            self._copy_repo("shared_files/", join(f, "shared_files/"))
            upload_command = ["shed_upload", "--tar_only"]
            upload_command.extend(self._shed_args())
            self._check_exit_code(upload_command)
            self._check_tar(
                f,
                "shed_upload.tar.gz",
                contains=[
                    "up_root/README.rst",
                    "up_root/cat.xml",
                    "shared_files/extra_test_data/extra_test_file.txt",
                ],
                not_contains=[],
            )

    def _assert_shed_diff(self, diff=1):
        shed_diff_command = ["shed_diff"]
        shed_diff_command.extend(self._shed_args())
        self._check_exit_code(shed_diff_command, exit_code=diff)

    def _verify_expansion(self, f, name=None):
        upload_command = ["shed_upload", "--tar_only"]
        upload_command.extend(self._shed_args())
        if name is not None:
            upload_command.append(join(f, name))
        self._check_exit_code(upload_command)
        self._check_tar(
            f,
            "shed_upload_cs_cat1.tar.gz",
            contains=["cat1.xml", "macros.xml", "test-data/1.bed"],
            not_contains=["cat2.xml"],
        )
        self._check_tar(
            f,
            "shed_upload_cs_cat2.tar.gz",
            contains=["cat2.xml", "macros.xml", "test-data/1.bed"],
            not_contains=["cat1.xml"],
        )

    def _verify_single_uploaded(self, f, download_args=None):
        self._verify_upload(f, ["cat.xml", "related_file", "test-data/1.bed"], download_args)

    def _verify_empty_repository(self, f, download_args=None):
        target = self._download_repo(f, download_args)
        assert len(os.listdir(target)) == 0

    def _verify_upload(self, f, download_files=None, download_args=None):
        download_files = download_files or []
        target = self._download_repo(f, download_args)
        for download_file in download_files:
            assert_exists(join(target, download_file))
        return target

    def _check_tar(self, f, tar_path, contains=None, not_contains=None):
        contains = contains or []
        not_contains = not_contains or []
        tar_path = join(f, tar_path)
        assert_exists(tar_path)
        target = self._untar(f, tar_path)
        for path in contains:
            assert_exists(join(target, path))
        for path in not_contains:
            assert not exists(join(target, path))
        return target

    def _download_repo(self, f, download_args=None):
        download_command = ["shed_download"]
        if download_args:
            download_command.extend(download_args)
        download_command.extend(self._shed_args(read_only=True))
        self._check_exit_code(download_command)
        download = join(f, "shed_download.tar.gz")
        assert_exists(download)
        return self._untar(f, "shed_download.tar.gz", tarbomb=False)

    def _untar(self, f, path, tarbomb=True):
        target = join(f, "download")
        if exists(target):
            shutil.rmtree(target)
        os.makedirs(target)
        try:
            tar = tarfile.open(path, "r:gz")
        except tarfile.ReadError as e:
            # Fixed in later version of Python, see
            # http://bugs.python.org/issue6123
            assert unicodify(e) == "empty header", e
            return target  # note contained no files!
        tar.extractall(path=target)
        tar.close()
        tar = tarfile.open(path, "r:gz")
        for tar_info in tar.getmembers():
            # These entries cause problems with TS.
            assert tar_info.name != "."
            assert tar_info.name != ""
        tar.close()
        if not tarbomb:
            return os.path.join(target, os.listdir(target)[0])
        else:
            return target


def update_package_1(f):
    """Update tool dependencies file for package_1."""
    changed_xml = join(TEST_REPOS_DIR, "package_1_changed", "tool_dependencies.xml")
    shutil.copyfile(changed_xml, join(f, "tool_dependencies.xml"))
