import os
from tempfile import mkstemp
import tarfile
import yaml
try:
    from bioblend import toolshed
except ImportError:
    toolshed = None


# Planemo generated or consumed files that do not need to be uploaded to the
# tool shed.
PLANEMO_FILES = [
    "tool_test_output.html",
    ".travis",
    ".travis.yml",
    ".shed.yml"
]
SHED_SHORT_NAMES = {
    "toolshed": "https://toolshed.g2.bx.psu.edu/",
    "testtoolshed": "https://testtoolshed.g2.bx.psu.edu/",
    "local": "http://localhost:9009/"
}


def shed_repo_config(path):
    shed_yaml_path = os.path.join(path, ".shed.yml")
    if os.path.exists(shed_yaml_path):
        with open(shed_yaml_path, "r") as f:
            return yaml.load(f)
    else:
        return {}


def tool_shed_client(ctx, **kwds):
    shed_target = kwds.get("shed_target")
    global_config = ctx.global_config
    if global_config and "sheds" in global_config:
        sheds_config = global_config["sheds"]
        shed_config = sheds_config.get(shed_target, {})
    else:
        shed_config = {}

    def prop(key):
        return kwds.get("shed_%s" % key, None) or shed_config.get(key, None)

    url = _tool_shed_url(kwds)
    key = prop("key")
    email = prop("email")
    password = prop("password")
    tsi = toolshed.ToolShedInstance(
        url=url,
        key=key,
        email=email,
        password=password
    )
    return tsi


def find_repository_id(ctx, tsi, path, **kwds):
    repo_config = shed_repo_config(path)
    owner = kwds.get("owner", None) or repo_config.get("owner", None)
    name = kwds.get("name", None) or repo_config.get("name", None)
    if owner is None:
        owner = ctx.global_config.get("shed_username", None)
    if name is None:
        name = os.path.basename(os.path.abspath(path))
    repos = tsi.repositories.get_repositories()

    def matches(r):
        return r["owner"] == owner and r["name"] == name

    matching_repos = filter(matches, repos)
    if not matching_repos:
        message = "Failed to find repository for owner/name %s/%s"
        raise Exception(message % (owner, name))
    repo_id = matching_repos[0]["id"]
    return repo_id


def build_tarball(tool_path):
    """Build a tool-shed tar ball for the specified path, caller is
    responsible for deleting this file.
    """
    fd, temp_path = mkstemp()
    try:
        with tarfile.open(temp_path, "w:gz") as tar:
            for name in os.listdir(tool_path):
                path = os.path.join(tool_path, name)
                tar.add(path, name, recursive=True, exclude=_tar_excludes)
    finally:
        os.close(fd)
    return temp_path


def _tool_shed_url(kwds):
    url = kwds.get("shed_target")
    if url in SHED_SHORT_NAMES:
        url = SHED_SHORT_NAMES[url]
    return url


def _tar_excludes(path):
    name = os.path.basename(path)
    if name.startswith(".git"):
        return True
    elif name in PLANEMO_FILES:
        return True
    return False
