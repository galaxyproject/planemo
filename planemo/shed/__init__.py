from collections import namedtuple
import contextlib
import copy
import fnmatch
import hashlib
import json
import os
import re
import shutil
import sys
import tarfile
from tempfile import (
    mkstemp,
)

from six import iteritems
import yaml

from planemo.io import (
    error,
    shell,
    info,
    can_write_to_path,
    temp_directory,
)
from planemo import git
from planemo import glob
from planemo.tools import load_tool_elements_from_path
from planemo import templates

from .interface import (
    username,
    tool_shed_instance,
    find_repository,
    api_exception_to_message,
    find_category_ids,
    download_tar,
)
from .diff import diff_and_remove

SHED_CONFIG_NAME = '.shed.yml'
REPO_DEPENDENCIES_CONFIG_NAME = "repository_dependencies.xml"
TOOL_DEPENDENCIES_CONFIG_NAME = "tool_dependencies.xml"

NO_REPOSITORIES_MESSAGE = ("Could not find any .shed.yml files or a --name to "
                           "describe the target repository.")
NAME_INVALID_MESSAGE = ("Cannot use --name argument when multiple directories "
                        "in target contain .shed.yml files.")
NAME_REQUIRED_MESSAGE = ("No repository name discovered but oneis required.")
CONFLICTING_NAMES_MESSAGE = ("The supplied name argument --name conflicts "
                             "with value discovered in .shed.yml.")
PARSING_PROBLEM = ("Problem parsing file .shed.yml in directory %s, skipping "
                   "repository. Message: [%s].")
AUTO_REPO_CONFLICT_MESSAGE = ("Cannot specify both auto_tool_repositories and "
                              "repositories in .shed.yml at this time.")
AUTO_NAME_CONFLICT_MESSAGE = ("Cannot specify both auto_tool_repositories and "
                              "in .shed.yml and --name on the command-line.")
REALIZAION_PROBLEMS_MESSAGE = ("Problem encountered executing action for one or more "
                               "repositories.")
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
REPO_TYPE_UNRESTRICTED = "unrestricted"
REPO_TYPE_TOOL_DEP = "tool_dependency_definition"
REPO_TYPE_SUITE = "repository_suite_definition"

# TODO: sync this with tool shed impl someday
VALID_REPOSITORYNAME_RE = re.compile("^[a-z0-9\_]+$")
VALID_PUBLICNAME_RE = re.compile("^[a-z0-9\-]+$")

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

    create_failed = _create_shed_config(ctx, shed_config_path, **kwds)
    if create_failed:
        return 1

    repo_dependencies_path = os.path.join(path, REPO_DEPENDENCIES_CONFIG_NAME)
    from_workflow = kwds.get("from_workflow", None)

    if from_workflow:
        workflow_name = os.path.basename(from_workflow)
        workflow_target = os.path.join(path, workflow_name)
        if not os.path.exists(workflow_target):
            shutil.copyfile(from_workflow, workflow_target)

        if not can_write_to_path(repo_dependencies_path, **kwds):
            return 1

        repo_pairs = _parse_repos_from_workflow(from_workflow)
        repository_dependencies = RepositoryDependencies()
        repository_dependencies.repo_pairs = repo_pairs
        repository_dependencies.write_to_path(repo_dependencies_path)

    return 0


def upload_repository(ctx, realized_repository, **kwds):
    """Upload a tool directory as a tarball to a tool shed.
    """
    path = realized_repository.path
    tar_path = kwds.get("tar", None)
    if not tar_path:
        tar_path = build_tarball(path, **kwds)
    if kwds["tar_only"]:
        name = realized_repository.pattern_to_file_name("shed_upload.tar.gz")
        shell("cp '%s' '%s'" % (tar_path, name))
        return 0
    tsi = tool_shed_client(ctx, **kwds)
    update_kwds = {}
    _update_commit_message(ctx, realized_repository, update_kwds, **kwds)

    repo_id = realized_repository.find_repository_id(ctx, tsi)
    if repo_id is None and kwds["force_repository_creation"]:
        repo_id = realized_repository.create(ctx, tsi)
    # failing to create the repo, give up
    if repo_id is None:
        return -1

    if kwds["check_diff"]:
        is_diff = diff_repo(ctx, realized_repository, **kwds)
        if not is_diff:
            name = realized_repository.name
            info("Repositry [%s] not different, skipping upload." % name)
            return 0

    # TODO: support updating repo information if it changes in the config file
    try:
        tsi.repositories.update_repository(repo_id, tar_path, **update_kwds)
    except Exception as e:
        message = api_exception_to_message(e)
        error("Could not update %s" % realized_repository.name)
        error(message)
        return -1
    info("Repository %s updated successfully." % realized_repository.name)
    return 0


def _update_commit_message(ctx, realized_repository, update_kwds, **kwds):
    message = kwds.get("message", None)
    git_rev = realized_repository.git_rev(ctx)
    git_repo = realized_repository.git_repo(ctx)
    if message is None:
        message = "planemo upload"
        if git_repo:
            message += " for repository %s" % git_repo
        if git_rev:
            message += " commit %s" % git_rev
    update_kwds["commit_message"] = message


def diff_repo(ctx, realized_repository, **kwds):
    with temp_directory("tool_shed_diff_") as working:
        return _diff_in(ctx, working, realized_repository, **kwds)


def _diff_in(ctx, working, realized_repository, **kwds):
    path = realized_repository.path
    shed_target_source = kwds.get("shed_target_source", None)

    label_a = "_%s_" % (shed_target_source if shed_target_source else "local")
    shed_target = kwds.get("shed_target", "B")
    if "/" in shed_target:
        shed_target = "custom_shed"
    label_b = "_%s_" % shed_target

    mine = os.path.join(working, label_a)
    other = os.path.join(working, label_b)

    tsi = tool_shed_client(ctx, read_only=True, **kwds)
    download_tarball(
        ctx,
        tsi,
        realized_repository,
        destination=other,
        clean=True,
        **kwds
    )
    if shed_target_source:
        new_kwds = kwds.copy()
        new_kwds["shed_target"] = shed_target_source
        tsi = tool_shed_client(ctx, read_only=True, **new_kwds)
        download_tarball(
            ctx,
            tsi,
            realized_repository,
            destination=mine,
            clean=True,
            **new_kwds
        )
    else:
        tar_path = build_tarball(path)
        cmd_template = 'mkdir "%s"; tar -xzf "%s" -C "%s"; rm -rf %s'
        shell(cmd_template % (mine, tar_path, mine, tar_path))

    output = kwds.get("output", None)
    raw = kwds.get("raw", False)
    is_diff = 0
    if not raw:
        if output:
            with open(output, "w") as f:
                is_diff = diff_and_remove(working, label_a, label_b, f)
        else:
            is_diff = diff_and_remove(working, label_a, label_b, sys.stdout)

    cmd = 'cd "%s"; diff -r %s %s' % (working, label_a, label_b)
    if output:
        cmd += ">> '%s'" % output
    exit = shell(cmd) or is_diff
    return exit


def shed_repo_config(path, name=None):
    shed_yaml_path = os.path.join(path, SHED_CONFIG_NAME)
    config = {}
    if os.path.exists(shed_yaml_path):
        with open(shed_yaml_path, "r") as f:
            config = yaml.load(f)

    if config is None:  # yaml may yield None
        config = {}
    _expand_raw_config(config, path, name=name)
    return config


def tool_shed_client(ctx=None, **kwds):
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
    return tool_shed_instance(url, key, email, password)


def find_repository_id(ctx, tsi, path, **kwds):
    repo_config = kwds.get("config", None)
    if repo_config is None:
        name = kwds.get("name", None)
        repo_config = shed_repo_config(path, name=name)
    name = repo_config["name"]
    find_kwds = kwds.copy()
    if "name" in find_kwds:
        del find_kwds["name"]
    return _find_repository_id(ctx, tsi, name, repo_config, **find_kwds)


def _find_repository_id(ctx, tsi, name, repo_config, **kwds):
    global_config = getattr(ctx, "global_config", {})
    owner = kwds.get("owner", None) or repo_config.get("owner", None)
    if owner is None:
        owner = global_config.get("shed_username", None)
    if owner is None:
        owner = username(tsi)

    matching_repository = find_repository(tsi, owner, name)
    if matching_repository is None:
        if not kwds.get("allow_none", False):
            message = "Failed to find repository for owner/name %s/%s"
            raise Exception(message % (owner, name))
        else:
            return None
    else:
        repo_id = matching_repository["id"]
        return repo_id


def _expand_raw_config(config, path, name=None):
    name_input = name
    if "name" not in config:
        config["name"] = name
    if config["name"] is None:
        config["name"] = path_to_repo_name(path)

    default_include = config.get("include", ["**"])
    repos = config.get("repositories", None)
    auto_tool_repos = config.get("auto_tool_repositories", False)
    suite_config = config.get("suite", False)

    if repos and auto_tool_repos:
        raise Exception(AUTO_REPO_CONFLICT_MESSAGE)
    if auto_tool_repos and name_input:
        raise Exception(AUTO_NAME_CONFLICT_MESSAGE)
    if auto_tool_repos:
        repos = _build_auto_tool_repos(path, config, auto_tool_repos)
    if suite_config:
        if repos is None:
            repos = {}
        _build_suite_repo(config, repos, suite_config)
    # If repositories aren't defined, just define a single
    # one based on calculated name and including everything
    # by default.
    if repos is None:
        repos = {
            config["name"]: {
                "include": default_include
            }
        }
    config["repositories"] = repos


def _build_auto_tool_repos(path, config, auto_tool_repos):
    default_include = config.get("include", ["**"])
    tool_els = list(load_tool_elements_from_path(path, recursive=True))
    paths = list(map(lambda pair: pair[0], tool_els))
    excludes = _shed_config_excludes(config)

    def _build_repository(tool_path, tool_el):
        tool_id = tool_el.getroot().get("id")
        tool_name = tool_el.getroot().get("name")
        template_vars = dict(
            tool_id=tool_id,
            tool_name=tool_name,
        )
        other_paths = paths[:]
        other_paths.remove(tool_path)
        tool_excludes = excludes + list(other_paths)
        repo_dict = {
            "include": default_include,
            "exclude": tool_excludes,
        }
        for key in ["name", "description", "long_description"]:
            template_key = "%s_template" % key
            template = auto_tool_repos.get(template_key, None)
            if template:
                value = templates.render(template, **template_vars)
                repo_dict[key] = value
        return repo_dict

    repos = {}
    for tool_path, tool_el in tool_els:
        repository_config = _build_repository(tool_path, tool_el)
        repository_name = repository_config["name"]
        repos[repository_name] = repository_config
    return repos


def _build_suite_repo(config, repos, suite_config):
    if not isinstance(suite_config, dict):
        suite_config = {}

    name = suite_config.get("name", None)
    if name is None:
        raise Exception("suite_configitories required name key.")
    description = suite_config.get("description", "")
    long_description = suite_config.get("long_description", None)
    owner = config["owner"]

    repo_pairs = map(lambda name: (owner, name), repos.keys())
    extra_repos = suite_config.get("include_repositories", {})
    extra_pairs = map(lambda item: (item["owner"], item["name"]), extra_repos)

    repository_dependencies = RepositoryDependencies()
    repository_dependencies.description = description
    repository_dependencies.repo_pairs = list(repo_pairs) + list(extra_pairs)

    repo = {
        "_files": {
            REPO_DEPENDENCIES_CONFIG_NAME: str(repository_dependencies)
        },
        "include": [],
        "name": name,
        "description": description,
    }
    if long_description:
        repo["long_description"] = long_description
    repos[name] = repo


def create_repository_for(ctx, tsi, name, repo_config):
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


def download_tarball(ctx, tsi, realized_repository, **kwds):
    repo_id = realized_repository.find_repository_id(ctx, tsi)
    if repo_id is None:
        message = "Unable to find repository id, cannot download."
        error(message)
        raise Exception(message)
    destination_pattern = kwds.get('destination', 'shed_download.tar.gz')
    destination = realized_repository.pattern_to_file_name(destination_pattern)
    to_directory = not destination.endswith("gz")
    download_tar(tsi, repo_id, destination, to_directory=to_directory)
    if to_directory:
        clean = kwds.get("clean", False)
        if clean:
            archival_file = os.path.join(destination, ".hg_archival.txt")
            if os.path.exists(archival_file):
                os.remove(archival_file)


def build_tarball(realized_path, **kwds):
    """Build a tool-shed tar ball for the specified path, caller is
    responsible for deleting this file.
    """

    # Simplest solution to sorting the files is to use a list,
    files = []
    for dirpath, dirnames, filenames in os.walk(realized_path):
        for f in filenames:
            files.append(os.path.join(dirpath, f))
    files.sort()

    fd, temp_path = mkstemp()
    try:
        tar = tarfile.open(temp_path, "w:gz", dereference=True)
        try:
            for raw in files:
                name = os.path.relpath(raw, realized_path)
                tar.add(os.path.join(realized_path, name), arcname=name)
        finally:
            tar.close()
    finally:
        os.close(fd)
    return temp_path


def for_each_repository(function, path, **kwds):
    ret_codes = []
    with _path_on_disk(path) as raw_path:
        try:
            for realized_repository in _realize_effective_repositories(
                raw_path, **kwds
            ):
                ret_codes.append(
                    function(realized_repository)
                )
        except RealizationException:
            error(REALIZAION_PROBLEMS_MESSAGE)
            return 254

    # "Good" returns are Nones, everything else is a -1 and should be
    # passed upwards.
    return 0 if all((not x) for x in ret_codes) else -1


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


def _realize_effective_repositories(path, **kwds):
    """ Expands folders in a source code repository into tool shed
    repositories.

    Each folder may have nested repositories and each folder may corresponding
    to many repositories (for instance if a folder has n tools in the source
    code repository but are published to the tool shed as one repository per
    tool).
    """
    raw_repo_objects = _find_raw_repositories(path, **kwds)
    failed = False
    with temp_directory() as base_dir:
        for raw_repo_object in raw_repo_objects:
            if isinstance(raw_repo_object, Exception):
                _handle_realization_error(raw_repo_object, **kwds)
                failed = True
                continue

            realized_repos = raw_repo_object.realizations(
                base_dir,
                kwds.get("fail_on_missing", True)
            )
            for realized_repo in realized_repos:
                if isinstance(realized_repo, Exception):
                    _handle_realization_error(realized_repo, **kwds)
                    failed = True
                    continue
                yield realized_repo
    if failed:
        raise RealizationException()


def _create_shed_config(ctx, path, **kwds):
    name = kwds.get("name", None) or path_to_repo_name(os.path.dirname(path))
    name_invalid = validate_repo_name(name)
    if name_invalid:
        error(name_invalid)
        return 1

    owner = kwds.get("owner", None)
    if owner is None:
        owner = ctx.global_config.get("shed_username", None)
    owner_invalid = validate_repo_owner(owner)
    if owner_invalid:
        error(owner_invalid)
        return 1
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


@contextlib.contextmanager
def _path_on_disk(path):
    git_path = None
    if path.startswith("git:"):
        git_path = path
    elif path.startswith("git+"):
        git_path = path[len("git+"):]
    if git_path is None:
        yield path
    else:
        with temp_directory() as git_repo:
            # TODO: pass ctx down through
            git.clone(None, git_path, git_repo)
            yield git_repo


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
        shed_file_dir = shed_file_dirs[0]
        try:
            config = shed_repo_config(shed_file_dir, name=name)
        except Exception as e:
            error_message = PARSING_PROBLEM % (shed_file_dir, e)
            return [RuntimeError(error_message)]
        config_name = config.get("name", None)

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
    multiple = len(raw_dirs) > 1
    name = kwds.get("name", None)

    # List of RawRepositoryDirectories or parsing failures if
    # fail_fast is not enabled.
    raw_repo_objects = []
    for raw_dir in raw_dirs:
        try:
            config = shed_repo_config(raw_dir, name=name)
        except Exception as e:
            error_message = PARSING_PROBLEM % (raw_dir, e)
            exception = RuntimeError(error_message)
            _handle_realization_error(exception, **kwds)
            raw_repo_objects.append(exception)
            continue
        raw_repo_object = RawRepositoryDirectory(raw_dir, config, multiple)
        raw_repo_objects.append(raw_repo_object)
    return raw_repo_objects


class RepositoryDependencies(object):
    """ Abstraction for shed repository_dependencies.xml files.
    """

    def __init__(self):
        self.description = ""
        self.repo_pairs = []

    def __str__(self):
        contents = '<repositories description="%s">' % self.description
        line_template = '  <repository owner="%s" name="%s" />\n'
        for (owner, name) in self.repo_pairs:
            contents += line_template % (owner, name)
        contents += "</repositories>"
        return contents

    def write_to_path(self, path):
        with open(path, "w") as f:
            f.write(str(self).encode("utf-8"))


class RawRepositoryDirectory(object):

    def __init__(self, path, config, multiple):
        self.path = path
        self.config = config
        self.name = config["name"]
        self.type = shed_repo_type(config, self.name)
        self.multiple = multiple  # operation over many repos?

    def _hash(self, name):
        return hashlib.md5(name.encode('utf-8')).hexdigest()

    def realizations(self, parent_directory, fail_on_missing=True):
        names = self._repo_names()

        for name in names:
            directory = os.path.join(parent_directory, self._hash(name), name)
            multiple = self.multiple or len(names) > 1
            if not os.path.exists(directory):
                os.makedirs(directory)
            yield self._realize_to(directory, name, multiple, fail_on_missing)

    def _realize_to(self, directory, name, multiple, fail_on_missing):
        ignore_list = []
        config = self._realize_config(name)
        excludes = _shed_config_excludes(config)
        for exclude in excludes:
            ignore_list.extend(_glob(self.path, exclude))

        realized_files = self._realized_files(name)
        missing = realized_files.include_failures
        if missing and fail_on_missing:
            msg = "Failed to include files for %s" % missing
            return RuntimeError(msg)

        for realized_file in realized_files.files:
            relative_dest = realized_file.relative_dest
            implicit_ignore = self._implicit_ignores(relative_dest)
            explicit_ignore = (realized_file.absolute_src in ignore_list)
            if implicit_ignore or explicit_ignore:
                continue
            realized_file.realize_to(directory)

        for (name, contents) in iteritems(config.get("_files", {})):
            path = os.path.join(directory, name)
            with open(path, "w") as f:
                f.write(contents)

        return RealizedRepositry(
            realized_path=directory,
            real_path=self.path,
            config=config,
            multiple=multiple,
            missing=missing,
        )

    def _repo_names(self):
        return self.config.get("repositories").keys()

    def _realized_files(self, name):
        config = self._realize_config(name)
        realized_files = []
        missing = []
        for include in config["include"]:
            included = RealizedFile.realized_files_for(self.path, include)
            if not included:
                missing.append(include)
            else:
                realized_files.extend(included)
        return RealizedFiles(realized_files, missing)

    def _realize_config(self, name):
        config = copy.deepcopy(self.config)
        config["name"] = name
        repo_config = config.get("repositories", {}).get(name, {})
        config.update(repo_config)
        if "repositories" in config:
            del config["repositories"]
        return config

    def _implicit_ignores(self, relative_path):
        # Filter out "unwanted files" :) like READMEs for special
        # repository types.
        if self.type == REPO_TYPE_TOOL_DEP:
            if relative_path != TOOL_DEPENDENCIES_CONFIG_NAME:
                return True

        if self.type == REPO_TYPE_SUITE:
            if relative_path != REPO_DEPENDENCIES_CONFIG_NAME:
                return True

        name = os.path.basename(relative_path)
        if relative_path.startswith(".git"):
            return True
        elif name in PLANEMO_FILES:
            return True
        return False

RealizedFiles = namedtuple("RealizedFiles", ["files", "include_failures"])


class RealizedFile(object):

    def __init__(self, src_root, src, dest, dest_is_file, strip_components):
        self.src_root = src_root
        self.src = src
        self.dest = dest
        self.dest_is_file = dest_is_file
        self.strip_components = strip_components

    @property
    def relative_dest(self):
        if self.dest_is_file:
            destination = self.dest
        else:
            destination = os.path.join(self.dest, self.stripped_source)
        return os.path.relpath(destination)

    @property
    def stripped_source(self):
        return "/".join(self.src.split("/")[self.strip_components:])

    @property
    def absolute_src(self):
        return os.path.abspath(os.path.join(self.src_root, self.src))

    def realize_to(self, directory):
        source_path = self.absolute_src
        if os.path.islink(source_path):
            source_path = os.path.realpath(source_path)
        relative_dest = self.relative_dest
        target_path = os.path.join(directory, relative_dest)
        target_exists = os.path.exists(target_path)
        if not target_exists:
            target_dir = os.path.dirname(target_path)
            if not os.path.exists(target_dir):
                os.makedirs(target_dir)
            if os.path.isdir(source_path):
                os.makedirs(target_path)
            else:
                os.symlink(source_path, target_path)

    @staticmethod
    def realized_files_for(path, include_info):
        if not isinstance(include_info, dict):
            include_info = {"source": include_info}
        source = include_info.get("source")
        destination = include_info.get("destination", None)
        strip_components = include_info.get("strip_components", 0)
        if destination is None:
            destination = "./"
            destination_specified = False
        else:
            destination_specified = True
        abs_source = os.path.join(path, source)
        dest_is_file = destination_specified and os.path.isfile(abs_source)
        realized_files = []
        for globbed_file in _glob(path, source):
            src = os.path.relpath(globbed_file, path)
            realized_files.append(
                RealizedFile(path, src, destination, dest_is_file, strip_components)
            )
        return realized_files

    def __str__(self):
        return "RealizedFile[src={},dest={},dest_file={}]".format(
            self.src, self.dest, self.dest_is_file
        )


class RealizedRepositry(object):

    def __init__(self, realized_path, real_path, config, multiple, missing):
        self.path = realized_path
        self.real_path = real_path
        self.config = config
        self.name = config["name"]
        self.multiple = multiple
        self.missing = missing

    def git_rev(self, ctx):
        return git.rev_if_git(ctx, self.real_path)

    def git_repo(self, ctx):
        return self.config.get("remote_repository_url", None)

    def pattern_to_file_name(self, pattern):
        if not self.multiple:
            return pattern

        name = self.config["name"]
        suffix = "_%s" % name.replace("-", "_")

        if "." not in pattern:
            return pattern + suffix
        else:
            parts = pattern.split(".", 1)
            return parts[0] + suffix + "." + parts[1]

    def find_repository_id(self, ctx, tsi):
        try:
            repo_id = _find_repository_id(
                ctx,
                tsi,
                name=self.name,
                repo_config=self.config,
                allow_none=True,
            )
            return repo_id
        except Exception as e:
            message = api_exception_to_message(e)
            error("Could not update %s" % self.name)
            error(message)
        return None

    def create(self, ctx, tsi):
        """Wrapper for creating the endpoint if it doesn't exist
        """
        try:
            repo = create_repository_for(
                ctx,
                tsi,
                self.name,
                self.config,
            )
            return repo['id']
        # Have to catch missing snyopsis/bioblend exceptions
        except Exception as e:
            # TODO: galaxyproject/bioblend#126
            try:
                upstream_error = json.loads(e.read())
                error(upstream_error['err_msg'])
            except Exception:
                error(str(e))
            return None


def _glob(path, pattern):
    pattern = os.path.join(path, pattern)
    if os.path.isdir(pattern):
        pattern = "%s/**" % pattern
    return glob.glob(pattern)


def _shed_config_excludes(config):
    return config.get('ignore', []) + config.get('exclude', [])


def _handle_realization_error(exception, **kwds):
    fail_fast = kwds.get("fail_fast", False)
    if fail_fast:
        raise exception
    else:
        error(str(exception))


def validate_repo_name(name):
    def _build_error(descript):
        return "Repository name [%s] invalid. %s" % (name, descript)

    msg = None
    if len(name) < 2:
        msg = _build_error(
            "Repository names must be at least 2 characters in length."
        )
    if len(name) > 80:
        msg = _build_error(
            "Repository names cannot be more than 80 characters in length."
        )
    if not VALID_REPOSITORYNAME_RE.match(name):
        msg = _build_error(
            "Repository names must contain only lower-case letters, "
            "numbers and underscore."
        )
    return msg


def validate_repo_owner(owner):
    def _build_error(descript):
        return "Owner [%s] invalid. %s" % (owner, descript)
    msg = None
    if len(owner) < 3:
        msg = _build_error("Owner must be at least 3 characters in length")
    if len(owner) > 255:
        msg = _build_error(
            "Owner cannot be more than 255 characters in length"
        )
    if not(VALID_PUBLICNAME_RE.match(owner)):
        msg = _build_error(
            "Owner must contain only lower-case letters, numbers and '-'"
        )
    return msg


class RealizationException(Exception):
    """ This exception indicates there was a problem while
    realizing effective repositories for a shed command. As a
    precondition - the user has already been informed with error().
    """

__all__ = [
    'for_each_repository',
    'api_exception_to_message',
    'tool_shed_client',
    'diff',
    'download_tarball',
    'shed_init',
    'CURRENT_CATEGORIES',
    'path_to_repo_name',
    'REPO_TYPE_UNRESTRICTED',
    'REPO_TYPE_TOOL_DEP',
    'REPO_TYPE_SUITE',
]
