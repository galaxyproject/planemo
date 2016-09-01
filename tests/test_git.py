import contextlib
import os

from planemo import git

from .test_utils import io

EXPECTED_HELLO_REV = "1c36390f585f8baa953548c00fc18c58e32fcf8b"

COMMITTER_DATE = "GIT_COMMITTER_DATE='2000-01-01T00:00:00+0000'"
COMMITTER_NAME = "GIT_COMMITTER_NAME='a' GIT_COMMITTER_EMAIL='a@example.com'"
COMMIT = ("git commit --date='2000-01-01T00:00:00+0000' "
          "--author='a <a@example.com>' -m 'initial'")


def test_rev():
    with _git_directory() as t:
        rev = git.rev(None, t)
        assert rev == EXPECTED_HELLO_REV, rev


def test_rev_if_git():
    with io.temp_directory() as t:
        rev = git.rev_if_git(None, t)
        assert rev is None


def test_rev_dirty_if_git():
    with _git_directory() as t:
        io.write_file(os.path.join(t, "README"), "Hello World!")
        rev = git.rev_if_git(None, t)
        assert rev == EXPECTED_HELLO_REV + "-dirty", rev


@contextlib.contextmanager
def _git_directory():
    with io.temp_directory() as t:
        io.write_file(os.path.join(t, "README"), "Hello!")
        cmd = " && ".join([
            "cd '{0}'",
            "git init .",
            "git add README",
            "%s %s %s" % (
                COMMITTER_DATE,
                COMMITTER_NAME,
                COMMIT
            )
        ]).format(t)
        assert io.shell(cmd) == 0
        yield t
