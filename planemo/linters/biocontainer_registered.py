"""Ensure a best-practice BioContainer is registered for a tool."""

from typing import (
    List,
    NamedTuple,
    Optional,
    TYPE_CHECKING,
)

from galaxy.tool_util.deps.container_resolvers.mulled import targets_to_mulled_name
from galaxy.tool_util.lint import Linter

from planemo.conda import tool_source_conda_targets
from planemo.linters.util import (
    cached_lint_result,
    xml_node_from_toolsource,
)

if TYPE_CHECKING:
    from galaxy.tool_util.deps.conda_util import CondaTarget
    from galaxy.tool_util.lint import LintContext
    from galaxy.tool_util.parser.interface import ToolSource

MESSAGE_WARN_NO_REQUIREMENTS = "No valid package requirement tags found to infer BioContainer from."
MESSAGE_WARN_NO_CONTAINER = "Failed to find a BioContainer registered for these requirements."
MESSAGE_INFO_FOUND_BIOCONTAINER = "BioContainer best-practice container found [%s]."

lint_tool_types = ["*"]


class BiocontainerCheck(NamedTuple):
    targets: List["CondaTarget"]
    name: Optional[str]


def _biocontainer_check(tool_source: "ToolSource") -> BiocontainerCheck:
    def calculate() -> BiocontainerCheck:
        targets = tool_source_conda_targets(tool_source)
        name = mulled_container_name("biocontainers", targets) if targets else None
        return BiocontainerCheck(targets, name)

    return cached_lint_result(tool_source, "biocontainer_registered", calculate)


class BiocontainerMissing(Linter):
    @classmethod
    def lint(cls, tool_source: "ToolSource", lint_ctx: "LintContext"):
        check = _biocontainer_check(tool_source)
        if check.targets and not check.name:
            requirements_node = xml_node_from_toolsource(tool_source, "requirements")
            lint_ctx.warn(MESSAGE_WARN_NO_CONTAINER, linter=cls.name(), node=requirements_node)


class BiocontainerRequirementsMissing(Linter):
    @classmethod
    def lint(cls, tool_source: "ToolSource", lint_ctx: "LintContext"):
        if not _biocontainer_check(tool_source).targets:
            requirements_node = xml_node_from_toolsource(tool_source, "requirements")
            lint_ctx.warn(MESSAGE_WARN_NO_REQUIREMENTS, linter=cls.name(), node=requirements_node)


class BiocontainerValid(Linter):
    @classmethod
    def lint(cls, tool_source: "ToolSource", lint_ctx: "LintContext"):
        check = _biocontainer_check(tool_source)
        if check.name:
            requirements_node = xml_node_from_toolsource(tool_source, "requirements")
            lint_ctx.info(
                MESSAGE_INFO_FOUND_BIOCONTAINER % check.name,
                linter=cls.name(),
                node=requirements_node,
            )


def mulled_container_name(namespace: str, targets: List["CondaTarget"]) -> Optional[str]:
    name = targets_to_mulled_name(targets=targets, hash_func="v2", namespace=namespace)
    if name:
        return f"quay.io/{namespace}/{name}"
    return None
