"""Utilities for interacting with Github."""
from __future__ import absolute_import

import os

from galaxy.tool_util.deps.commands import which

from planemo import git
from planemo.io import (
    communicate,
    IS_OS_X,
    untar_to,
)

try:
    import github
    has_github_lib = True
except ImportError:
    github = None
    has_github_lib = False

GH_VERSION = "1.0.0"

NO_GITHUB_DEP_ERROR = ("Cannot use github functionality - "
                       "PyGithub library not available.")
FAILED_TO_DOWNLOAD_GH = "No gh executable available and it could not be installed."


def get_github_config(ctx, allow_anonymous=False):
    """Return a :class:`planemo.github_util.GithubConfig` for given configuration."""
    global_github_config = _get_raw_github_config(ctx)
    return GithubConfig(global_github_config, allow_anonymous=allow_anonymous)


def clone_fork_branch(ctx, target, path, **kwds):
    """Clone, fork, and branch a repository ahead of building a pull request."""
    git.checkout(
        ctx,
        target,
        path,
        branch=kwds.get("branch", None),
        remote="origin",
        from_branch="master"
    )
    if kwds.get("fork"):
        try:
            fork(ctx, path, **kwds)
        except Exception:
            pass


def fork(ctx, path, **kwds):
    """Fork the target repository using ``hub``."""
    gh_path = ensure_gh(ctx, **kwds)
    gh_env = get_gh_env(ctx, path, **kwds)
    cmd = [gh_path, "repo", "fork"]
    communicate(cmd, env=gh_env)


def get_or_create_repository(ctx, organization, name, **kwds):
    pass


def create_release(ctx, repository, version, **kwds):
    repository = get_or_create_repository()


def pull_request(ctx, path, message=None, **kwds):
    """Create a pull request against the origin of the path using ``gh``."""
    gh_path = ensure_gh(ctx, **kwds)
    gh_env = get_gh_env(ctx, path, **kwds)
    cmd = [gh_path, "pr", "create"]
    if message is None:
        cmd.append('--fill')
    else:
        lines = message.splitlines()
        cmd.extend(['--title', lines[0]])
        if len(lines) > 1:
            cmd.extend(["--body", "\n".join(lines[1:])])
    communicate(cmd, env=gh_env)


def get_gh_env(ctx, path, **kwds):
    """Return a environment dictionary to run hub with given user and repository target."""
    env = git.git_env_for(path).copy()
    github_config = _get_raw_github_config(ctx)
    if github_config is not None:
        if "access_token" in github_config:
            env["GITHUB_TOKEN"] = github_config["access_token"]

    return env


def ensure_gh(ctx, **kwds):
    """Ensure ``hub`` is on the system ``PATH``.

    This method will ensure ``hub`` is installed if it isn't available.

    For more information on ``hub`` checkout ...
    """
    gh_path = which("gh")
    if not gh_path:
        planemo_gh_path = os.path.join(ctx.workspace, "gh")
        if not os.path.exists(planemo_gh_path):
            _try_download_hub(planemo_gh_path)

        if not os.path.exists(planemo_gh_path):
            raise Exception(FAILED_TO_DOWNLOAD_GH)

        gh_path = planemo_gh_path
    return gh_path


def _try_download_gh(planemo_gh_path):
    link = _gh_link()
    # Strip URL base and .tgz at the end.
    basename = link.split("/")[-1].rsplit(".", 1)[0]
    untar_to(link, tar_args=['-zxvf', '-', "%s/bin/gh" % basename], path=planemo_gh_path)
    communicate(["chmod", "+x", planemo_gh_path])


def _get_raw_github_config(ctx):
    """Return a :class:`planemo.github_util.GithubConfig` for given configuration."""
    if "github" not in ctx.global_config:
        if "GITHUB_TOKEN" in os.environ:
            return {
                "access_token": os.environ["GITHUB_TOKEN"],
            }
    if "github" not in ctx.global_config:
        raise Exception("github account not found in planemo config and GITHUB_TOKEN environment variables unset")
    return ctx.global_config["github"]


class GithubConfig(object):
    """Abstraction around a Github account.

    Required to use ``github`` module methods that require authorization.
    """

    def __init__(self, config, allow_anonymous=False):
        if not has_github_lib:
            raise Exception(NO_GITHUB_DEP_ERROR)
        if "access_token" not in config:
            if not allow_anonymous:
                raise Exception("github authentication unavailable")
            github_object = github.Github()
        else:
            github_object = github.Github(config["access_token"])
        self._github = github_object


def _gh_link():
    if IS_OS_X:
        template_link = "https://github.com/cli/cli/releases/download/v%s/gh_%s_macOS_amd64.tgz"
    else:
        template_link = "https://github.com/cli/cli/releases/download/v%s/gh_%s_linux_amd64.tgz"
    return template_link % (GH_VERSION, GH_VERSION)


def publish_as_gist_file(ctx, path, name="index"):
    """Publish a gist.

    More information on gists at http://gist.github.com/.
    """
    github_config = get_github_config(ctx, allow_anonymous=False)
    user = github_config._github.get_user()
    with open(path, "r") as fh:
        content = fh.read()
    content_file = github.InputFileContent(content)
    gist = user.create_gist(False, {name: content_file})
    return gist.files[name].raw_url


def get_repository_object(ctx, name):
    github_object = get_github_config(ctx, allow_anonymous=True)
    return github_object._github.get_repo(name)


__all__ = (
    "clone_fork_branch",
    "ensure_gh",
    "fork",
    "get_github_config",
    "get_gh_env",
    "publish_as_gist_file",
)
