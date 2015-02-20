import os
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

    warned = False
    for required_key in ["owner", "name"]:
        if required_key not in shed_contents:
            lint_ctx.warn(".shed.yml did not contain key [%s]" % required_key)
            warned = True

    if not warned:
        lint_ctx.info(".shed.yml found and appears to be valid YAML.")
