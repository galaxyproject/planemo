"""
"""

try:
    import github
    has_github_lib = True
except ImportError:
    github = None
    has_github_lib = False

NO_GITHUB_DEP_ERROR = ("Cannot use github functionality - "
                       "PyGithub library not available.")


def get_github_config(ctx):
    if "github" not in ctx.global_config:
        return None
    global_github_config = ctx.global_config["github"]
    return GithubConfig(global_github_config)


class GithubConfig(object):

    def __init__(self, config):
        if not has_github_lib:
            raise Exception(NO_GITHUB_DEP_ERROR)
        self._github = github.Github(config["username"], config["password"])


def publish_as_gist_file(ctx, path, name="index"):
    github_config = get_github_config(ctx)
    user = github_config._github.get_user()
    content = open(path, "r").read()
    content_file = github.InputFileContent(content)
    gist = user.create_gist(False, {name: content_file})
    return gist.files[name].raw_url
