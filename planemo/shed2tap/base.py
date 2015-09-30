from __future__ import print_function

from xml.etree import ElementTree

from six.moves.urllib.request import urlretrieve
from six.moves.urllib.error import URLError
from six import string_types

import os
import tarfile
import zipfile
from ftplib import all_errors as FTPErrors  # tuple of exceptions

TOOLSHED_MAP = {
    "toolshed": "https://toolshed.g2.bx.psu.edu",
    "testtoolshed": "https://testtoolshed.g2.bx.psu.edu",
}


class Dependencies(object):
    """ Base class for parsing Tool Shed dependency files.
    """

    def __init__(
        self,
        dependencies_file,
        repo=None,
        package_factory=None,
    ):
        if package_factory is None:
            package_factory = BasePackage
        self.repo = repo
        self.root = ElementTree.parse(dependencies_file).getroot()
        packages = []
        dependencies = []
        package_els = self.root.findall("package")
        assert package_els is not None
        for package_el in package_els:
            install_els = package_el.findall("install")
            readme_els = package_el.findall("readme")
            if len(readme_els) > 0:
                readme = readme_els[0].text
            else:
                readme = None
            assert len(install_els) in (0, 1)
            if len(install_els) == 1:
                install_el = install_els[0]
                package = package_factory(
                    self,
                    package_el,
                    install_el,
                    readme=readme
                )
                packages.append(package)
            else:
                repository_el = package_el.find("repository")
                if repository_el is None:
                    message = "no repository in package el for %s" % repo
                    raise AssertionError(message)
                dependency = Dependency(self, package_el, repository_el)
                dependencies.append(dependency)

        self.packages = packages
        self.dependencies = dependencies

    def single_package(self):
        return len(self.packages) == 1

    def __repr__(self):
        return "Dependencies[for_repo=%s]" % self.repo


class Repo(object):

    def __init__(self, **kwds):
        for key, value in kwds.iteritems():
            setattr(self, key, value)

    def recipe_base_name(self):
        owner = self.owner.replace("-", "")
        name = self.name
        name = name.replace("_", "").replace("-", "")
        base = "%s_%s" % (owner, name)
        return base

    @staticmethod
    def from_xml(elem):
        tool_shed_url = elem.attrib.get("toolshed", None)
        if tool_shed_url and ("testtoolshed" in tool_shed_url):
            prefix = "testtoolshed"
        else:
            prefix = "toolshed"
        prior = elem.attrib.get("prior_installation_required", False)
        return Repo(
            prefix=prefix,
            name=elem.attrib["name"],
            owner=elem.attrib["owner"],
            tool_shed_url=tool_shed_url,
            changeset_revision=elem.attrib.get("changeset_revision", None),
            prior_installation_required=prior,
        )

    @staticmethod
    def from_api(prefix, repo_json):
        return Repo(
            prefix=prefix,
            name=repo_json["name"],
            owner=repo_json["owner"],
            tool_shed_url=TOOLSHED_MAP[prefix],
        )

    def get_file(self, path):
        try:
            url_template = "%s/repos/%s/%s/raw-file/tip/%s"
            url = url_template % (
                self.tool_shed_url,
                self.owner,
                self.name,
                path
            )
            path, headers = urlretrieve(url)
            return path
        except Exception as e:
            print(e)
            return None

    def __repr__(self):
        return "Repository[name=%s,owner=%s]" % (self.name, self.owner)


class Dependency(object):

    def __init__(self, dependencies, package_el, repository_el):
        self.dependencies = dependencies
        self.package_el = package_el
        self.repository_el = repository_el
        self.repo = Repo.from_xml(repository_el)

    def __repr__(self):
        temp = "Dependency[package_name=%s,version=%s,dependent_package=%s]"
        return temp % (
            self.package_el.attrib["name"],
            self.package_el.attrib["version"],
            self.repository_el.attrib["name"]
        )


class BasePackage(object):

    def __init__(self, dependencies, package_el, install_el, readme):
        self.dependencies = dependencies
        self.package_el = package_el
        self.install_el = install_el
        self.readme = readme
        self.all_actions = self.get_all_actions()
        self.no_arch_option = self.has_no_achitecture_install()

    def get_all_actions(self):
        action_or_group = self.install_el[0]
        parsed_actions = []
        if action_or_group.tag == "actions":
            parsed_actions.append(self.parse_actions(action_or_group))
        elif action_or_group.tag == "actions_group":
            actions_els = action_or_group.findall("actions")
            assert actions_els is not None
            for actions in actions_els:
                parsed_actions.append(self.parse_actions(actions))
            action_els = action_or_group.findall("action")
            assert action_els is not None
            for action in action_els:
                for parsed_a in parsed_actions:
                    parsed_a.actions.append(self.parse_action(action))
        return parsed_actions

    def has_no_achitecture_install(self):
        all_actions = self.all_actions
        if len(all_actions) < 2:
            return False
        else:
            last_action = all_actions[-1]
            return (not last_action.architecture) and (not last_action.os)

    def has_explicit_set_environments(self):
        all_actions = self.all_actions
        for actions in all_actions:
            for action in actions.actions:
                if action.explicit_variables:
                    return True
        return False

    def has_multiple_set_environments(self):
        all_actions = self.all_actions
        for actions in all_actions:
            count = 0
            for action in actions.actions:
                if action.explicit_variables:
                    count += 1
            if count > 1:
                return True
        return False

    def parse_actions(self, actions):
        os = actions.attrib.get("os", None)
        architecture = actions.get("architecture", None)
        action_els = actions.findall("action")
        assert action_els is not None
        parsed_actions = map(self.parse_action, action_els)
        action_packages = []
        for package in actions.findall("package"):
            action_packages.append(self.parse_action_package(package))
        return Actions(parsed_actions, os, architecture, action_packages)

    def parse_action_package(self, elem):
        name = elem.attrib["name"]
        version = elem.attrib["version"]
        repo = Repo.from_xml(elem.find("repository"))
        return ActionPackage(name, version, repo)

    def parse_action(self, action):
        return BaseAction.from_elem(action, package=self)

    def __repr__(self):
        actions = self.all_actions
        parts = (
            self.package_el.attrib["name"],
            self.package_el.attrib["version"],
            self.dependencies,
            actions
        )
        template = "Install[name=%s,version=%s,dependencies=%s,actions=%s]"
        return template % parts


class Actions(object):

    def __init__(
        self,
        actions,
        os=None,
        architecture=None,
        action_packages=[]
    ):
        self.os = os
        self.architecture = architecture
        self.actions = actions or []
        self.action_packages = action_packages

    def first_download(self):
        for action in self.actions:
            if action.action_type in ["download_by_url", "download_file"]:
                return action
        return None

    def downloads(self):
        actions = []
        for action in self.actions:
            if action.action_type in ["download_by_url", "download_file"]:
                actions.append(action)
        return actions

    def __repr__(self):
        platform = ""
        if self.os or self.architecture:
            platform = "os=%s,arch=%s," % (self.os, self.architecture)
        return "Actions[%s%s]" % (platform, map(str, self.actions))

    def _indent_extend(self, target, new_entries, indent="    "):
        for line in new_entries:
            target.append(indent + line)

    def to_bash(self):
        # Use self.os.title() to match "Linux" or "Darwin" in bash where case matters:
        if self.os and self.architecture:
            condition = '("%s" == `uname`) && ("%s" == `arch`)' % (self.os.title(), self.architecture)
        elif self.os:
            condition = '"%s" == `uname`' % self.os.title()
        elif self.architecture:
            condition = '"%s" == `arch`' % self.architecture
        else:
            condition = None

        install_cmds = []
        env_cmds = []

        if condition:
            # Conditional actions block
            install_cmds = [
                '#' + '-' * 60,
                'if [[ $specifc_action_done == 0 && %s ]]' % condition,
                'then',
                '    echo "Platform-specific action for os=%s, arch=%s"' % (self.os, self.architecture)]
            env_cmds = install_cmds[:]
            # TODO - Refactor block indentation?
            for action in self.actions:
                i_cmds, e_cmds = action.to_bash()
                self._indent_extend(install_cmds, i_cmds)
                self._indent_extend(env_cmds, e_cmds)
            # If we run the action, do not want to run any later actions!
            install_cmds.extend(['    specifc_action_done=1', 'fi'])
            env_cmds.extend(['    specifc_action_done=1', 'fi'])
        else:
            # Non-specific default action...
            install_cmds = [
                '#' + '-' * 60,
                'if [[ $specifc_action_done == 0 ]]',
                'then',
                '    echo "Non-platform-specific actions"']
            env_cmds = install_cmds[:]
            for action in self.actions:
                i_cmds, e_cmds = action.to_bash()
                self._indent_extend(install_cmds, i_cmds)
                self._indent_extend(env_cmds, e_cmds)
            install_cmds.append('fi')
            env_cmds.append('fi')
        return install_cmds, env_cmds


class ActionPackage(object):

    def __init__(self, name, version, repo):
        self.name = name
        self.version = version
        self.repo = repo


class BaseAction(object):

    def __repr__(self):
        return "Action[type=%s]" % self.action_type

    def same_as(self, other):
        if self._keys != other._keys:
            return False
        else:
            for key in self._keys:
                if getattr(self, key) != getattr(other, key):
                    return False

            return True

    def parse_action_repo(self, elem):
        repo_elem = elem.find("repository")
        repo = Repo.from_xml(repo_elem)
        self.repo = repo

    def parse_package_elems(self, elem):
        package_els = elem.findall("package")
        packages = []
        assert package_els is not None
        for package_el in package_els:
            packages.append(package_el.text)
        self.packages = packages

    @classmethod
    def from_elem(cls, elem, package):
        type = elem.attrib["type"]
        action_class = actions_by_type[type]
        return action_class(elem)

    def to_bash(self):
        """Return lists of bash shell commands to execute this action.

        This method is be implemented by each sub-class, and will
        return two list of strings (for ``dep_install.sh`` and
        ``env.sh`` respectively).
        """
        raise NotImplementedError("No to_bash defined for %r" % self)


def _tar_folders(filename):
    archive = tarfile.open(filename, "r", errorlevel=0)
    folders = set()
    for i in archive.getmembers():
        if i.isdir():
            folders.add(i.name.rstrip("/"))
        else:
            folders.add(os.path.split(i.name)[0])
    return folders


def _zip_folders(filename):
    archive = zipfile.ZipFile(filename, "r")
    return set(i.filename.rstrip("/") for i in archive.infolist() if i.filename.endswith("/"))


def _common_prefix(folders):
    common_prefix = ""
    if len(folders) == 1:
        common_prefix = list(folders)[0]
    else:
        common_prefix = os.path.commonprefix(folders)
        assert not os.path.isabs(common_prefix), folders
    return common_prefix


def _cache_download(url, filename):
    """Returns local path to cached copy of URL using given filename."""
    try:
        cache = os.environ["DOWNLOAD_CACHE"]
    except KeyError:
        # TODO - expose this as a command line option
        raise ValueError("Dependencies cache location $DOWNLOAD_CACHE not set.")

    local = os.path.join(cache, filename)

    if not os.path.isfile(local):
        # Must download it...
        try:
            import sys  # TODO - log this nicely...
            sys.stderr.write("Downloading %s\n" % url)
            urlretrieve(url, local)
        except URLError:
            # Most likely server is down, could be bad URL in XML action:
            raise RuntimeError("Unable to download %s" % url)
        except FTPErrors:
            # Most likely server is down, could be bad URL in XML action:
            raise RuntimeError("Unable to download %s" % url)

    return local


def _determine_compressed_file_folder(url, downloaded_filename):
    """Determine how to decompress the file & its directory structure.

    Returns a list of shell commands. Consider this example where the
    folder to change to cannot be guessed from the tar-ball filename:

        $ curl -o "ncbi-blast-2.2.30+-ia32-linux.tar.gz" \
        "ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.30/ncbi-blast-2.2.30+-ia32-linux.tar.gz"
        $ tar -zxvf ncbi-blast-2.2.30+-ia32-linux.tar.gz
        $ cd ncbi-blast-2.2.30+

    Here it would return:

        ['tar -zxvf ncbi-blast-2.2.30+-ia32-linux.tar.gz', 'cd ncbi-blast-2.2.30+']

    If not cached, this function will download the file to the
    $DOWNLOAD_CACHE folder, and then open it / decompress it in
    order to find common folder prefix used.  This will also verify
    how to decompress the file.
    """
    answer = []
    local = _cache_download(url, downloaded_filename)

    if tarfile.is_tarfile(local):
        folders = _tar_folders(local)
        if downloaded_filename.endswith(".tar.gz") or downloaded_filename.endswith(".tgz"):
            answer.append('tar -zxvf %s' % downloaded_filename)
        elif downloaded_filename.endswith(".tar.bz2"):
            answer.append('tar -jxvf %s' % downloaded_filename)
        elif downloaded_filename.endswith(".tar"):
            answer.extend('tar -xvf %s' % downloaded_filename)
        else:
            # Quite possibly this file doesn't need decompressing,
            # but until we've tested lots of real world tool_dependencies.xml
            # files I'd like to check these cases to confirm this.
            raise NotImplementedError("How to decompress tar file %s?" % downloaded_filename)
    elif zipfile.is_zipfile(local):
        if local.endswith('.jar'):
            # Do not decompress!
            return answer
        folders = _zip_folders(local)
        answer.append('unzip %s' % downloaded_filename)
    elif downloaded_filename.endswith(".dmg"):
        # Do not decompress!
        return answer
    else:
        # No compression? Leave as it is?
        raise NotImplementedError("What kind of compression is %s using?" % local)

    common_prefix = _common_prefix(folders)
    if common_prefix:
        answer.append('cd "%s"' % common_prefix)

    return answer


def _commands_to_download_and_extract(url, target_filename=None):
    # TODO - Include checksum validation here?
    if target_filename:
        downloaded_filename = target_filename
    else:
        downloaded_filename = os.path.split(url)[-1]
        if "?" in downloaded_filename:
            downloaded_filename = downloaded_filename[:downloaded_filename.index("?")]
        if "#" in downloaded_filename:
            downloaded_filename = downloaded_filename[:downloaded_filename.index("#")]

    # Curl is present on Mac OS X, can we assume it will be on Linux?
    # Cannot assume that wget will be on Mac OS X.
    answer = [
        'if [[ -f "%s" ]]' % downloaded_filename,
        'then',
        '    echo "Reusing existing %s"' % downloaded_filename,
        'elif [[ -f "$DOWNLOAD_CACHE/%s" ]]' % downloaded_filename,
        'then',
        '    echo "Reusing cached %s"' % downloaded_filename,
        '    ln -s "$DOWNLOAD_CACHE/%s" "%s"' % (downloaded_filename, downloaded_filename),
        'else',
        '    echo "Downloading %s"' % downloaded_filename,
        '    curl -L -o "$DOWNLOAD_CACHE/%s" "%s"' % (downloaded_filename, url),
        '    ln -s "$DOWNLOAD_CACHE/%s" "%s"' % (downloaded_filename, downloaded_filename),
        'fi',
        ]
    answer.extend(_determine_compressed_file_folder(url, downloaded_filename))
    return answer, []


class DownloadByUrlAction(BaseAction):
    action_type = "download_by_url"
    _keys = ["url"]

    def __init__(self, elem):
        self.url = elem.text.strip()
        assert self.url

    def to_bash(self):
        # See class DownloadByUrl in Galaxy,
        # lib/tool_shed/galaxy_install/tool_dependencies/recipe/step_handler.py
        # Do we need to worry about target_filename here?
        return _commands_to_download_and_extract(self.url)


class DownloadFileAction(BaseAction):
    action_type = "download_file"
    _keys = ["url", "extract"]

    def __init__(self, elem):
        self.url = elem.text.strip()
        self.extract = asbool(elem.attrib.get("extract", False))

    def to_bash(self):
        if self.extract:
            return _commands_to_download_and_extract(self.url)
        else:
            return ['wget %s' % self.url], []


class DownloadBinary(BaseAction):
    action_type = "download_binary"
    _keys = ["url_template", "target_directory"]

    def __init__(self, elem):
        self.url_template = elem.text
        assert self.url_template
        self.target_directory = elem.get('target_directory', None)


class ShellCommandAction(BaseAction):
    action_type = "shell_command"
    _keys = ["command"]

    def __init__(self, elem):
        self.command = elem.text

    def to_bash(self):
        return [self.command], []


class TemplateShellCommandAction(BaseAction):
    action_type = "template_command"
    _keys = ["language", "command"]

    def __init__(self, elem):
        self.command = elem.text
        self.language = elem.get('language', 'cheetah').lower()
        assert self.command
        assert self.language == "cheetah"


class MoveFileAction(BaseAction):
    action_type = "move_file"
    _keys = ["move_file"]

    def __init__(self, elem):
        self.source = elem.find("source").text
        self.destination = elem.find("destination").text

    def to_bash(self):
        return ["mv %s %s" % (self.source, self.destination)], []


class MoveDirectoryFilesAction(BaseAction):
    action_type = "move_directory_files"
    _keys = ["source_directory", "destination_directory"]

    def __init__(self, elem):
        source_directory = elem.find("source_directory").text
        destination_directory = elem.find("destination_directory").text
        self.source_directory = source_directory
        self.destination_directory = destination_directory

    def to_bash(self):
        return ["mv %s/* %s/" % (self.source_directory, self.destination_directory)], []


class SetEnvironmentAction(BaseAction):
    action_type = "set_environment"
    _keys = ["variables"]

    def __init__(self, elem):
        variables = []
        var_els = elem.findall("environment_variable")
        assert var_els is not None
        for ev_elem in var_els:
            var = SetVariable(ev_elem)
            variables.append(var)
        self.variables = variables
        assert self.variables

    def to_bash(self):
        answer = []
        for var in self.variables:
            # Expand $INSTALL_DIR here?
            if var.action == "set_to":
                answer.append('export %s=%s' % (var.name, var.raw_value))
            elif var.action == "prepend_to":
                answer.append('export %s=%s:$%s' % (var.name, var.raw_value, var.name))
            elif var.action == "append_to":
                answer.append('export %s=$%s:%s' % (var.name, var.name, var.raw_value))
            else:
                raise ValueError("Undefined environment variable action %r" % var.action)
        return answer, answer  # Actions needed in env.sh here!


class ChmodAction(BaseAction):
    action_type = "chmod"
    _keys = ["mods"]

    def __init__(self, elem):
        mods = []
        file_els = elem.findall("file")
        assert file_els is not None
        for mod_elem in file_els:
            mod = {}
            mod["mode"] = mod_elem.attrib["mode"]
            mod["target"] = mod_elem.text
            mods.append(mod)
        self.mods = mods
        assert self.mods

    def to_bash(self):
        return ["chmod %s %s" % (m["mode"], m["target"]) for m in self.mods], []


class MakeInstallAction(BaseAction):
    action_type = "make_install"
    _keys = []

    def __init__(self, elem):
        pass

    def to_bash(self):
        return ["make install"], []


class AutoconfAction(BaseAction):
    action_type = "autoconf"
    _keys = ["options"]

    def __init__(self, elem):
        self.options = elem.text

    def to_bash(self):
        if self.options:
            raise NotImplementedError("Options with action autoconf not implemented yet.")
        return ['./configure', 'make', 'make install'], []


class ChangeDirectoryAction(BaseAction):
    action_type = "change_directory"
    _keys = ["directory"]

    def __init__(self, elem):
        self.directory = elem.text
        assert self.directory

    def to_bash(self):
        return ["cd %s" % self.directory], []


class MakeDirectoryAction(BaseAction):
    action_type = "make_directory"
    _keys = ["directory"]

    def __init__(self, elem):
        self.directory = elem.text

    def to_bash(self):
        return ["mkdir -p %s" % self.directory], []


class SetupPerlEnvironmentAction(BaseAction):
    action_type = "setup_perl_environment"
    _keys = ["repo", "packages"]

    def __init__(self, elem):
        self.parse_action_repo(elem)
        self.parse_package_elems(elem)


class SetupRubyEnvironmentAction(BaseAction):
    action_type = "setup_ruby_environment"
    _keys = ["repo", "packages"]

    def __init__(self, elem):
        self.parse_action_repo(elem)
        self.parse_package_elems(elem)


class SetupPythonEnvironmentAction(BaseAction):
    action_type = "setup_python_environment"
    _keys = ["repo", "packages"]

    def __init__(self, elem):
        self.parse_action_repo(elem)
        self.parse_package_elems(elem)


class SetupVirtualenvAction(BaseAction):
    action_type = "setup_virtualenv"
    _keys = ["use_requirements_file", "python", "requirements"]

    def __init__(self, elem):
        use_reqs = elem.attrib.get("use_requirements_file", "True")
        self.use_requirements_file = asbool(use_reqs)
        self.python = elem.get('python', 'python')
        self.requirements = elem.text or 'requirements.txt'


class SetupREnvironmentAction(BaseAction):
    action_type = "setup_r_environment"
    _keys = ["repo", "packages"]

    def __init__(self, elem):
        self.parse_action_repo(elem)
        self.parse_package_elems(elem)


class SetEnvironmentForInstallAction(BaseAction):
    action_type = "set_environment_for_install"

    def __init__(self, elem):
        pass

    def to_bash(self):
        # TODO - How could we resolve/check the dependencies?
        return ['echo "WARNING: Assuming packages already installed!"'], []


class SetVariable(object):

    def __init__(self, elem):
        self.action = elem.attrib["action"]
        self.name = elem.attrib["name"]
        self.raw_value = elem.text


truthy = frozenset(['true', 'yes', 'on', 'y', 't', '1'])
falsy = frozenset(['false', 'no', 'off', 'n', 'f', '0'])


def asbool(obj):
    if isinstance(obj, string_types):
        obj = obj.strip().lower()
        if obj in truthy:
            return True
        elif obj in falsy:
            return False
        else:
            raise ValueError("String is not true/false: %r" % obj)
    return bool(obj)


action_classes = [
    DownloadByUrlAction,
    DownloadFileAction,
    DownloadBinary,
    ShellCommandAction,
    TemplateShellCommandAction,
    MoveFileAction,
    MoveDirectoryFilesAction,
    SetEnvironmentAction,
    ChmodAction,
    MakeInstallAction,
    AutoconfAction,
    ChangeDirectoryAction,
    MakeDirectoryAction,
    SetupPerlEnvironmentAction,
    SetupRubyEnvironmentAction,
    SetupPythonEnvironmentAction,
    SetupVirtualenvAction,
    SetupREnvironmentAction,
    SetEnvironmentForInstallAction,
]

actions_by_type = dict(map(lambda c: (c.action_type, c), action_classes))
