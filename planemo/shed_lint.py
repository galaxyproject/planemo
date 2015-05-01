import os
import yaml
from galaxy.tools.lint import LintContext
from galaxy.tools.linters.help import rst_invalid
from planemo.lint import lint_xsd
from planemo.shed import (
    path_to_repo_name,
    REPO_TYPE_UNRESTRICTED,
    REPO_TYPE_TOOL_DEP,
    REPO_TYPE_SUITE,
    CURRENT_CATEGORIES,
    validate_repo_owner,
    validate_repo_name,
)
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


VALID_REPOSITORY_TYPES = [
    REPO_TYPE_UNRESTRICTED,
    REPO_TYPE_TOOL_DEP,
    REPO_TYPE_SUITE,
]


def lint_repository(ctx, realized_repository, **kwds):
    # TODO: this really needs to start working with realized path.
    path = realized_repository.real_path
    info("Linting repository %s" % path)
    lint_args = build_lint_args(ctx, **kwds)
    lint_ctx = LintContext(lint_args["level"])
    lint_ctx.lint(
        "lint_expansion",
        lint_expansion,
        realized_repository,
    )

    lint_ctx.lint(
        "lint_tool_dependencies",
        lint_tool_dependencies,
        path,
    )
    lint_ctx.lint(
        "lint_repository_dependencies",
        lint_repository_dependencies,
        path,
    )
    lint_ctx.lint(
        "lint_shed_yaml",
        lint_shed_yaml,
        path,
    )
    lint_ctx.lint(
        "lint_readme",
        lint_readme,
        path,
    )
    if kwds["tools"]:
        for (tool_path, tool_xml) in yield_tool_xmls(ctx, path,
                                                     recursive=True):
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


def lint_expansion(realized_repository, lint_ctx):
    missing = realized_repository.missing
    if missing:
        msg = "Failed to expand inclusions %s" % missing
        lint_ctx.warn(msg)
    else:
        lint_ctx.info("Included files all found.")


def lint_readme(path, lint_ctx):
    readme_rst = os.path.join(path, "README.rst")
    readme = os.path.join(path, "README")
    readme_txt = os.path.join(path, "README.txt")

    readme_found = False
    for readme in [readme_rst, readme, readme_txt]:
        if os.path.exists(readme):
            readme_found = readme

    readme_md = os.path.join(path, "README.md")
    if not readme_found and os.path.exists(readme_md):
        lint_ctx.warn("Tool Shed doesn't render markdown, "
                      "README.md is invalid readme.")
        return

    if not readme_found:
        # TODO: filter on TYPE and make this a warning if
        # unrestricted repository - need to update iuc standards
        # first though.
        lint_ctx.info("No README found skipping.")
        return

    if readme_found.endswith(".rst"):
        readme_text = open(readme_found, "r").read()
        invalid_rst = rst_invalid(readme_text)
        if invalid_rst:
            template = "Invalid restructured text found in README [%s]."
            msg = template % invalid_rst
            lint_ctx.warn(msg)
            return
        lint_ctx.info("README found containing valid reStructuredText.")
    else:
        lint_ctx.info("README found containing plain text.")


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
        return
    lint_ctx.info(".shed.yml found and appears to be valid YAML.")
    _lint_shed_contents(lint_ctx, path, shed_contents)


def _lint_shed_contents(lint_ctx, path, shed_contents):
    name = shed_contents.get("name", None)
    effective_name = name or path_to_repo_name(path)

    def _lint_if_present(key, func, *args):
        value = shed_contents.get(key, None)
        if value is not None:
            msg = func(value, *args)
            if msg:
                lint_ctx.warn(msg)

    _lint_if_present("owner", validate_repo_owner)
    _lint_if_present("name", validate_repo_name)
    _lint_if_present("type", _validate_repo_type, effective_name)
    _lint_if_present("categories", _validate_categories)


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


def _validate_categories(categories):
    msg = None
    if len(categories) == 0:
        msg = "Repository should specify one or more categories."
    else:
        for category in categories:
            unknown_categories = []
            if category not in CURRENT_CATEGORIES:
                unknown_categories.append(category)
            if unknown_categories:
                msg = "Categories [%s] unknown." % unknown_categories

    return msg
