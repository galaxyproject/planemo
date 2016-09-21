"""Utilities for interacting with Github."""
from __future__ import absolute_import

import os

from galaxy.tools.deps.commands import which

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

HUB_VERSION = "2.2.8"

NO_GITHUB_DEP_ERROR = ("Cannot use github functionality - "
                       "PyGithub library not available.")
FAILED_TO_DOWNLOAD_HUB = "No hub executable available and it could not be installed."


def get_github_config(ctx):
    """Return a :class:`planemo.github_util.GithubConfig` for given configuration."""
    global_github_config = _get_raw_github_config(ctx)
    return None if global_github_config is None else GithubConfig(global_github_config)


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
        fork(ctx, path, **kwds)


def fork(ctx, path, **kwds):
    """Fork the target repository using ``hub``."""
    hub_path = ensure_hub(ctx, **kwds)
    hub_env = get_hub_env(ctx, path, **kwds)
    cmd = [hub_path, "fork"]
    communicate(cmd, env=hub_env)


def pull_request(ctx, path, message=None, **kwds):
    """Create a pull request against the origin of the path using ``hub``."""
    hub_path = ensure_hub(ctx, **kwds)
    hub_env = get_hub_env(ctx, path, **kwds)
    cmd = [hub_path, "pull-request"]
    if message is not None:
        cmd.extend(["-m", message])
    communicate(cmd, env=hub_env)


def get_hub_env(ctx, path, **kwds):
    """Return a environment dictionary to run hub with given user and repository target."""
    env = git.git_env_for(path).copy()
    github_config = _get_raw_github_config(ctx)
    if github_config is not None:
        if "username" in github_config:
            env["GITHUB_USER"] = github_config["username"]
        if "password" in github_config:
            env["GITHUB_PASSWORD"] = github_config["password"]

    return env


def ensure_hub(ctx, **kwds):
    """Ensure ``hub`` is on the system ``PATH``.

    This method will ensure ``hub`` is installed if it isn't available.

    For more information on ``hub`` checkout ...
    """
    hub_path = which("hub")
    if not hub_path:
        planemo_hub_path = os.path.join(ctx.workspace, "hub")
        if not os.path.exists(planemo_hub_path):
            _try_download_hub(planemo_hub_path)

        if not os.path.exists(planemo_hub_path):
            raise Exception(FAILED_TO_DOWNLOAD_HUB)

        hub_path = planemo_hub_path
    return hub_path


def _try_download_hub(planemo_hub_path):
    link = _hub_link()
    # Strip URL base and .tgz at the end.
    basename = link.split("/")[-1].rsplit(".", 1)[0]
    untar_to(link, tar_args="-zxvf - %s/bin/hub -O > '%s'" % (basename, planemo_hub_path))
    communicate(["chmod", "+x", planemo_hub_path])


def _get_raw_github_config(ctx):
    """Return a :class:`planemo.github_util.GithubConfig` for given configuration."""
    if "github" not in ctx.global_config:
        return None
    return ctx.global_config["github"]


class GithubConfig(object):
    """Abstraction around a Github account.

    Required to use ``github`` module methods that require authorization.
    """

    def __init__(self, config):
        if not has_github_lib:
            raise Exception(NO_GITHUB_DEP_ERROR)
        self._github = github.Github(config["username"], config["password"])


def _hub_link():
    if IS_OS_X:
        template_link = "https://github.com/github/hub/releases/download/v%s/hub-darwin-amd64-%s.tgz"
    else:
        template_link = "https://github.com/github/hub/releases/download/v%s/hub-linux-amd64-%s.tgz"
    return template_link % (HUB_VERSION, HUB_VERSION)


def publish_as_gist_file(ctx, path, name="index"):
    """Publish a gist.

    More information on gists at http://gist.github.com/.
    """
    github_config = get_github_config(ctx)
    user = github_config._github.get_user()
    content = open(path, "r").read()
    content_file = github.InputFileContent(content)
    gist = user.create_gist(False, {name: content_file})
    return gist.files[name].raw_url


__all__ = [
    "clone_fork_branch",
    "ensure_hub",
    "fork",
    "get_github_config",
    "get_hub_env",
    "publish_as_gist_file",
]
