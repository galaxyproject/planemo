"""Ensure best-practice biocontainer registered for this tool."""

from typing import (
    List,
    Optional,
    TYPE_CHECKING,
)

from galaxy.tool_util.deps.container_resolvers.mulled import targets_to_mulled_name
from galaxy.tool_util.deps.mulled.mulled_build_tool import requirements_to_mulled_targets
from galaxy.tool_util.lint import Linter

from .util import xml_node_from_toolsource

if TYPE_CHECKING:
    from galaxy.tool_util.deps.conda_util import CondaTarget
    from galaxy.tool_util.lint import LintContext
    from galaxy.tool_util.parser.interface import ToolSource

MESSAGE_WARN_NO_REQUIREMENTS = "No valid package requirement tags found to infer BioContainer from."
MESSAGE_WARN_NO_CONTAINER = "Failed to find a BioContainer registered for these requirements."
MESSAGE_INFO_FOUND_BIOCONTAINER = "BioContainer best-practice container found [%s]."

lint_tool_types = ["*"]


class BiocontainerValid(Linter):
    @classmethod
    def lint(cls, tool_source: "ToolSource", lint_ctx: "LintContext"):
        requirements, *_ = tool_source.parse_requirements_and_containers()
        targets = requirements_to_mulled_targets(requirements)
        name = mulled_container_name("biocontainers", targets)
        if name:
            requirements_node = xml_node_from_toolsource(tool_source, "requirements")
            lint_ctx.info(MESSAGE_INFO_FOUND_BIOCONTAINER % name, linter=cls.name(), node=requirements_node)


class BiocontainerMissing(Linter):
    @classmethod
    def lint(cls, tool_source: "ToolSource", lint_ctx: "LintContext"):
        requirements, *_ = tool_source.parse_requirements_and_containers()
        targets = requirements_to_mulled_targets(requirements)
        name = mulled_container_name("biocontainers", targets)
        if not name:
            requirements_node = xml_node_from_toolsource(tool_source, "requirements")
            lint_ctx.warn(MESSAGE_WARN_NO_CONTAINER, linter=cls.name(), node=requirements_node)


def mulled_container_name(namespace: str, targets: List["CondaTarget"]) -> Optional[str]:
    name = targets_to_mulled_name(targets=targets, hash_func="v2", namespace=namespace)
    if name:
        return f"quay.io/{namespace}/{name}"
    else:
        return None
