import os
import re
import yaml
from galaxy.tools.lint import LintContext
from planemo.lint import lint_xsd
from planemo.tool_lint import (
    build_lint_args,
    yield_tool_xmls,
)
from planemo.xml import XSDS_PATH


from planemo.io import info
from planemo.io import error

from galaxy.tools.lint import lint_xml_with

TOOL_DEPENDENCIES_XSD = os.path.join(XSDS_PATH, "tool_dependencies.xsd")
REPO_DEPENDENCIES_XSD = os.path.join(XSDS_PATH, "repository_dependencies.xsd")

# TODO: sync this with tool shed impl someday
VALID_REPOSITORYNAME_RE = re.compile("^[a-z0-9\_]+$")
VALID_PUBLICNAME_RE = re.compile("^[a-z0-9\-]+$")

VALID_REPOSITORY_TYPES = [
    "unrestricted",
    "tool_dependency_definition",
    "repository_suite_definition",
]


def lint_repository(ctx, path, **kwds):
    info("Linting repository %s" % path)
    lint_args = build_lint_args(**kwds)
    lint_ctx = LintContext(lint_args["level"])
    lint_ctx.lint(
        "tool_dependencies",
        lint_tool_dependencies,
        path,
    )
    lint_ctx.lint(
        "repository_dependencies",
        lint_repository_dependencies,
        path,
    )
    lint_ctx.lint(
        "shed_yaml",
        lint_shed_yaml,
        path,
    )
    if kwds["tools"]:
        for (tool_path, tool_xml) in yield_tool_xmls(ctx, path):
            info("+Linting tool %s" % tool_path)
            lint_xml_with(
                lint_ctx,
                tool_xml,
                extra_modules=lint_args["extra_modules"]
            )
    failed = lint_ctx.failed(lint_args["fail_level"])
    if failed:
        error("Failed linting")
    return 1 if failed else 0


def lint_tool_dependencies(path, lint_ctx):
    tool_dependencies = os.path.join(path, "tool_dependencies.xml")
    if not os.path.exists(tool_dependencies):
        lint_ctx.info("No tool_dependencies.xml, skipping.")
        return
    lint_xsd(lint_ctx, TOOL_DEPENDENCIES_XSD, tool_dependencies)


def lint_repository_dependencies(path, lint_ctx):
    repo_dependencies = os.path.join(path, "repository_dependencies.xml")
    if not os.path.exists(repo_dependencies):
        lint_ctx.info("No repository_dependencies.xml, skipping.")
        return
    lint_xsd(lint_ctx, REPO_DEPENDENCIES_XSD, repo_dependencies)


def lint_shed_yaml(path, lint_ctx):
    shed_yaml = os.path.join(path, ".shed.yml")
    if not os.path.exists(shed_yaml):
        lint_ctx.info("No .shed.yml file found, skipping.")
        return
    try:
        shed_contents = yaml.load(open(shed_yaml, "r"))
    except Exception as e:
        lint_ctx.warn("Failed to parse .shed.yml file [%s]" % str(e))
    lint_ctx.info(".shed.yml found and appears to be valid YAML.")
    _lint_shed_contents(lint_ctx, path, shed_contents)


def _lint_shed_contents(lint_ctx, path, shed_contents):
    def _lint_if_present(key, func, *args):
        value = shed_contents.get(key, None)
        if value is not None:
            msg = func(value, *args)
            if msg:
                lint_ctx.warn(msg)

    _lint_if_present("owner", _validate_repo_owner)
    _lint_if_present("name", _validate_repo_name)
    effective_name = shed_contents.get("name", None) or os.path.basename(path)
    _lint_if_present("type", _validate_repo_type, effective_name)


def _validate_repo_type(repo_type, name):
    if repo_type not in VALID_REPOSITORY_TYPES:
        return "Invalid repository type specified [%s]" % repo_type

    is_dep = repo_type == "tool_dependency_definition"
    is_suite = repo_type == "repository_suite_definition"
    if is_dep and not name.startswith("package_"):
        return ("Tool dependency definition repositories should have names "
                "starting with package_")
    if is_suite and not name.startswith("suite_"):
        return ("Repository suite definition repositories should have names "
                "starting with suite_")
    if name.startswith("package_") or name.startswith("suite_"):
        if repo_type == "unrestricted":
            return ("Repository name indicated specialized repository type "
                    "but repository is listed as unrestricted.")


def _validate_repo_name(name):
    msg = None
    if len(name) < 2:
        msg = "Repository names must be at least 2 characters in length."
    if len(name) > 80:
        msg = "Repository names cannot be more than 80 characters in length."
    if not(VALID_REPOSITORYNAME_RE.match(name)):
        msg = ("Repository names must contain only lower-case letters, "
               "numbers and underscore.")
    return msg


def _validate_repo_owner(owner):
    msg = None
    if len(owner) < 3:
        msg = "Owner must be at least 3 characters in length"
    if len(owner) > 255:
        msg = "Owner cannot be more than 255 characters in length"
    if not(VALID_PUBLICNAME_RE.match(owner)):
        msg = "Owner must contain only lower-case letters, numbers and '-'"
    return msg
