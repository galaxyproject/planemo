import fnmatch
import glob
import hashlib
import json
import os
import tarfile
from tempfile import (
    mkstemp,
    mkdtemp,
)
import shutil

import yaml

try:
    from bioblend import toolshed
except ImportError:
    toolshed = None

from planemo.io import (
    error,
    untar_to,
    can_write_to_path,
)

SHED_CONFIG_NAME = '.shed.yml'
NO_REPOSITORIES_MESSAGE = ("Could not find any .shed.yml files or a --name to "
                           "describe the target repository.")
NAME_INVALID_MESSAGE = ("Cannot use --name argument when multiple directories "
                        "in target contain .shed.yml files.")
NAME_REQUIRED_MESSAGE = ("No repository name discovered but oneis required.")
CONFLICTING_NAMES_MESSAGE = ("The supplied name argument --name conflicts "
                             "with value discovered in .shed.yml.")
PARSING_PROBLEM = ("Problem parsing file .shed.yml in directory %s, skipping "
                   "repository. Message: [%s].")

# Planemo generated or consumed files that do not need to be uploaded to the
# tool shed.
PLANEMO_FILES = [
    "shed_upload.tar.gz",
    "tool_test_output.json",
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
REPOSITORY_DOWNLOAD_TEMPLATE = (
    "%srepository/download?repository_id=%s"
    "&changeset_revision=default&file_type=gz"
)
BIOBLEND_UNAVAILABLE = ("This functionality requires the bioblend library "
                        " which is unavailable, please install `pip install "
                        "bioblend`")

REPO_TYPE_UNRESTRICTED = "unrestricted"
REPO_TYPE_TOOL_DEP = "tool_dependency_definition"
REPO_TYPE_SUITE = "repository_suite_definition"

# Generate with python scripts/categories.py
CURRENT_CATEGORIES = [
    "Assembly",
    "ChIP-seq",
    "Combinatorial Selections",
    "Computational chemistry",
    "Convert Formats",
    "Data Managers",
    "Data Source",
    "Fasta Manipulation",
    "Fastq Manipulation",
    "Genome-Wide Association Study",
    "Genomic Interval Operations",
    "Graphics",
    "Imaging",
    "Metabolomics",
    "Metagenomics",
    "Micro-array Analysis",
    "Next Gen Mappers",
    "Ontology Manipulation",
    "Phylogenetics",
    "Proteomics",
    "RNA",
    "SAM",
    "Sequence Analysis",
    "Statistics",
    "Systems Biology",
    "Text Manipulation",
    "Tool Dependency Packages",
    "Tool Generators",
    "Transcriptomics",
    "Variant Analysis",
    "Visualization",
    "Web Services",
]


def shed_init(ctx, path, **kwds):
    if not os.path.exists(path):
        os.makedirs(path)
    shed_config_path = os.path.join(path, SHED_CONFIG_NAME)
    if not can_write_to_path(shed_config_path, **kwds):
        # .shed.yml exists and no --force sent.
        return 1

    _create_shed_config(ctx, shed_config_path, **kwds)
    repo_dependencies_path = os.path.join(path, "repository_dependencies.xml")
    from_workflow = kwds.get("from_workflow", None)

    if from_workflow:
        workflow_name = os.path.basename(from_workflow)
        workflow_target = os.path.join(path, workflow_name)
        if not os.path.exists(workflow_target):
            shutil.copyfile(from_workflow, workflow_target)

        if not can_write_to_path(repo_dependencies_path, **kwds):
            return 1

        repo_pairs = _parse_repos_from_workflow(from_workflow)
        contents = '<repositories description="">'
        line_template = '  <repository owner="%s" name="%s" />'
        for (owner, name) in repo_pairs:
            contents += line_template % (owner, name)
        contents += "</repositories>"
        with open(repo_dependencies_path, "w") as f:
            f.write(contents)

    return 0


def shed_repo_config(path):
    shed_yaml_path = os.path.join(path, SHED_CONFIG_NAME)
    if os.path.exists(shed_yaml_path):
        with open(shed_yaml_path, "r") as f:
            return yaml.load(f)
    else:
        return {}


def tool_shed_client(ctx=None, **kwds):
    if toolshed is None:
        raise Exception(BIOBLEND_UNAVAILABLE)
    read_only = kwds.get("read_only", False)
    shed_target = kwds.get("shed_target")
    global_config = getattr(ctx, "global_config", {})
    if global_config and "sheds" in global_config:
        sheds_config = global_config["sheds"]
        shed_config = sheds_config.get(shed_target, {})
    else:
        shed_config = {}

    def prop(key):
        return kwds.get("shed_%s" % key, None) or shed_config.get(key, None)

    url = _tool_shed_url(kwds)
    if read_only:
        key = None
        email = None
        password = None
    else:
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
    global_config = getattr(ctx, "global_config", {})
    repo_config = shed_repo_config(path)
    owner = kwds.get("owner", None) or repo_config.get("owner", None)
    name = kwds.get("name", None) or repo_config.get("name", None)
    if owner is None:
        owner = global_config.get("shed_username", None)
    if name is None:
        name = path_to_repo_name(path)
    repos = tsi.repositories.get_repositories()

    def matches(r):
        return r["owner"] == owner and r["name"] == name

    matching_repos = list(filter(matches, repos))
    if not matching_repos:
        if not kwds.get("allow_none", False):
            message = "Failed to find repository for owner/name %s/%s"
            raise Exception(message % (owner, name))
        else:
            return None
    else:
        repo_id = matching_repos[0]["id"]
        return repo_id


def create_repository(ctx, tsi, path, **kwds):
    repo_config = shed_repo_config(path)
    name = kwds.get("name", None) or repo_config.get("name", None)
    if name is None:
        name = path_to_repo_name(path)

    description = repo_config.get("description", None)
    long_description = repo_config.get("long_description", None)
    repo_type = shed_repo_type(repo_config, name)
    remote_repository_url = repo_config.get("remote_repository_url", None)
    homepage_url = repo_config.get("homepage_url", None)
    categories = repo_config.get("categories", [])
    category_ids = find_category_ids(tsi, categories)

    # description is required, as is name.
    if description is None:
        message = "description required for automatic creation of repositories"
        raise Exception(message)

    repo = tsi.repositories.create_repository(
        name=name,
        synopsis=description,
        description=long_description,
        type=repo_type,
        remote_repository_url=remote_repository_url,
        homepage_url=homepage_url,
        category_ids=category_ids
    )
    return repo


def find_category_ids(tsi, categories):
    """ Translate human readable category names into their associated IDs.
    """
    category_list = tsi.repositories.get_categories()

    category_ids = []
    for cat in categories:
        matching_cats = [x for x in category_list if x['name'] == cat]
        if not matching_cats:
            message = "Failed to find category %s" % cat
            raise Exception(message)
        category_ids.append(matching_cats[0]['id'])
    return category_ids


def download_tarball(ctx, tsi, path, **kwds):
    destination = kwds.get('destination', 'shed_download.tar.gz')
    repo_id = find_repository_id(ctx, tsi, path, **kwds)
    base_url = tsi.base_url
    if not base_url.endswith("/"):
        base_url += "/"
    download_url = REPOSITORY_DOWNLOAD_TEMPLATE % (base_url, repo_id)
    to_directory = not destination.endswith("gz")
    if to_directory:
        untar_args = "-xzf - -C %s --strip-components 1" % destination
    else:
        untar_args = None
    untar_to(download_url, destination, untar_args)
    if to_directory:
        clean = kwds.get("clean", False)
        if clean:
            archival_file = os.path.join(destination, ".hg_archival.txt")
            if os.path.exists(archival_file):
                os.remove(archival_file)


def build_tarball(tool_path, **kwds):
    """Build a tool-shed tar ball for the specified path, caller is
    responsible for deleting this file.
    """

    # Not really how realize_effective_repositories was meant to be used.
    # It should be pushed up a level into the thing that is uploading tar
    # balls to iterate over them - but placing it here for now because
    # it address some bugs.
    effective_repositories = realize_effective_repositories(tool_path, **kwds)
    for realized_repository in effective_repositories:
        fd, temp_path = mkstemp()
        try:
            tar = tarfile.open(temp_path, "w:gz")
            try:
                tar.add(realized_repository.path, arcname=".", recursive=True)
            finally:
                tar.close()
        finally:
            os.close(fd)
        return temp_path
    raise Exception("Problem not valid repositories found.")


def walk_repositories(path):
    """ Recurse through directories and find effective repositories. """
    for base_path, dirnames, filenames in os.walk(path):
        for filename in fnmatch.filter(filenames, '.shed.yml'):
            yield base_path


def for_each_repository(function, path):
    ret_codes = []
    for base_path in walk_repositories(path):
        ret_codes.append(
            function(base_path)
        )
    # "Good" returns are Nones, everything else is a -1 and should be
    # passed upwards.
    return 0 if all((not x) for x in ret_codes) else -1


def username(tsi):
    user = _user(tsi)
    return user["username"]


def _user(tsi):
    """ Fetch user information from the ToolShed API for given
    key.
    """
    # TODO: this should be done with an actual bioblend method,
    # see https://github.com/galaxyproject/bioblend/issues/130.
    response = tsi.make_get_request(tsi.url + "/users")
    return response.json()[0]


def path_to_repo_name(path):
    return os.path.basename(os.path.abspath(path))


def shed_repo_type(config, name):
    repo_type = config.get("type", None)
    if repo_type is None and name.startswith("package_"):
        repo_type = REPO_TYPE_TOOL_DEP
    elif repo_type is None and name.startswith("suite_"):
        repo_type = REPO_TYPE_SUITE
    elif repo_type is None:
        repo_type = REPO_TYPE_UNRESTRICTED
    return repo_type


def _tool_shed_url(kwds):
    url = kwds.get("shed_target")
    if url in SHED_SHORT_NAMES:
        url = SHED_SHORT_NAMES[url]
    return url


def realize_effective_repositories(path, **kwds):
    """ Expands folders in a source code repository into tool shed
    repositories.

    Each folder may have nested repositories and each folder may corresponding
    to many repositories (for instance if a folder has n tools in the source
    code repository but are published to the tool shed as one repository per
    tool).
    """
    raw_repo_objects = _find_raw_repositories(path, **kwds)
    temp_directory = mkdtemp()
    try:
        for raw_repo_object in raw_repo_objects:
            yield raw_repo_object.realize_to(temp_directory)
    finally:
        shutil.rmtree(temp_directory)


def _create_shed_config(ctx, path, **kwds):
    name = kwds.get("name", None) or path_to_repo_name(os.path.dirname(path))
    owner = kwds.get("owner", None)
    if owner is None:
        owner = ctx.global_config.get("shed_username", None)
    description = kwds.get("description", None) or name
    long_description = kwds.get("long_description", None)
    remote_repository_url = kwds.get("remote_repository_url", None)
    homepage_url = kwds.get("homepage_url", None)
    categories = kwds.get("category", [])
    config = dict(
        name=name,
        owner=owner,
        description=description,
        long_description=long_description,
        remote_repository_url=remote_repository_url,
        homepage_url=homepage_url,
        categories=categories,
    )
    # Remove empty entries...
    for k in list(config.keys()):
        if config[k] is None:
            del config[k]

    with open(path, "w") as f:
        yaml.dump(config, f)


def _parse_repos_from_workflow(path):
    with open(path, "r") as f:
        workflow_json = json.load(f)
    steps = workflow_json["steps"]
    tool_ids = set()
    for value in steps.values():
        step_type = value["type"]
        if step_type != "tool":
            continue
        tool_id = value["tool_id"]
        if "/repos/" in tool_id:
            tool_ids.add(tool_id)

    repo_pairs = set()
    for tool_id in tool_ids:
        tool_repo_info = tool_id.split("/repos/", 1)[1]
        tool_repo_parts = tool_repo_info.split("/")
        owner = tool_repo_parts[0]
        name = tool_repo_parts[1]
        repo_pairs.add((owner, name))

    return repo_pairs


def _find_raw_repositories(path, **kwds):
    name = kwds.get("name", None)
    recursive = kwds.get("recursive", False)

    shed_file_dirs = []
    if recursive:
        for base_path, dirnames, filenames in os.walk(path):
            for filename in fnmatch.filter(filenames, SHED_CONFIG_NAME):
                shed_file_dirs.append(base_path)
    elif os.path.exists(os.path.join(path, SHED_CONFIG_NAME)):
        shed_file_dirs.append(path)

    config_name = None
    if len(shed_file_dirs) == 1:
        config = shed_repo_config(shed_file_dirs[0])
        config_name = config.get("name", None)

    if len(shed_file_dirs) < 2 and config_name is None and name is None:
        name = path_to_repo_name(path)

    if len(shed_file_dirs) > 1 and name is not None:
        raise Exception(NAME_INVALID_MESSAGE)
    if config_name is not None and name is not None:
        if config_name != name:
            raise Exception(CONFLICTING_NAMES_MESSAGE)
    raw_dirs = shed_file_dirs or [path]
    kwds_copy = kwds.copy()
    kwds_copy["name"] = name
    return _build_raw_repo_objects(raw_dirs, **kwds_copy)


def _build_raw_repo_objects(raw_dirs, **kwds):
    """
    From specific directories with .shed.yml files or specified directly from
    the comman-line build abstract description of directories that should be
    expanded out into shed repositories.
    """
    name = kwds.get("name", None)
    skip_errors = kwds.get("skip_errors", False)

    raw_repo_objects = []
    for raw_dir in raw_dirs:
        try:
            config = shed_repo_config(raw_dir)
        except Exception as e:
            if skip_errors:
                error_message = PARSING_PROBLEM % (raw_dir, e)
                error(error_message)
            else:
                raise
        if name:
            config["name"] = name
        raw_repo_object = RawRepositoryDirectory(raw_dir, config)
        raw_repo_objects.append(raw_repo_object)
    return raw_repo_objects


class RawRepositoryDirectory(object):

    def __init__(self, path, config):
        self.path = path
        self.config = config
        self.name = config["name"]
        self.type = shed_repo_type(config, self.name)

    @property
    def _hash(self):
        return hashlib.md5(self.name.encode('utf-8')).hexdigest()

    def realize_to(self, parent_directory):
        directory = os.path.join(parent_directory, self._hash, self.name)
        if not os.path.exists(directory):
            os.makedirs(directory)

        ignore_list = []
        for shed_ignore in self.config.get('ignore', []):
            ignore_list.extend(glob.glob(os.path.join(self.path, shed_ignore)))

        for root, _, files in os.walk(self.path, followlinks=True):
            for name in files:
                full_path = os.path.join(root, name)
                relative_path = os.path.relpath(full_path, self.path)
                implicit_ignore = self._implicit_ignores(relative_path)
                explicit_ignore = (full_path in ignore_list)
                if implicit_ignore or explicit_ignore:
                    continue

                self._realize_file(relative_path, directory)
        return RealizedRepositry(directory, self.config)

    def _realize_file(self, relative_path, directory):
        source_path = os.path.join(self.path, relative_path)
        if os.path.islink(source_path):
            source_path = os.path.realpath(source_path)
        target_path = os.path.join(directory, relative_path)
        target_dir = os.path.dirname(target_path)
        if not os.path.exists(target_dir):
            os.makedirs(target_dir)
        os.symlink(source_path, target_path)

    def _implicit_ignores(self, relative_path):
        # Filter out "unwanted files" :) like READMEs for special
        # repository types.
        if self.type == REPO_TYPE_TOOL_DEP:
            if relative_path != "tool_dependencies.xml":
                return True

        if self.type == REPO_TYPE_SUITE:
            if relative_path != "repository_dependencies.xml":
                return True

        name = os.path.basename(relative_path)
        if relative_path.startswith(".git"):
            return True
        elif name in PLANEMO_FILES:
            return True
        return False


class RealizedRepositry(object):

    def __init__(self, path, config):
        self.path = path
        self.config = config
