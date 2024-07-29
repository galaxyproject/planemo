"""Ensure requirements are matched in best practice conda channels."""

from typing import (
    Generator,
    TYPE_CHECKING,
)

from galaxy.tool_util.deps.conda_util import requirement_to_conda_targets
from galaxy.tool_util.lint import Linter

from .util import xml_node_from_toolsource
from planemo.conda import (
    BEST_PRACTICE_CHANNELS,
    best_practice_search,
)

if TYPE_CHECKING:
    from galaxy.tool_util.lint import LintContext
    from galaxy.tool_util.parser.interface import ToolSource

lint_tool_types = ["*"]


class CondaRequirementValid(Linter):
    @classmethod
    def lint(cls, tool_source: "ToolSource", lint_ctx: "LintContext"):
        for conda_target in _requirements_conda_targets(tool_source):
            (best_hit, exact) = best_practice_search(conda_target)
            conda_target_str = conda_target.package
            if conda_target.version:
                conda_target_str += "@%s" % (conda_target.version)
            if best_hit and exact:
                message = f"Requirement [{conda_target_str}] matches target in best practice Conda channel [{best_hit.get('channel')}]."
                requirements_node = xml_node_from_toolsource(tool_source, "requirements")
                lint_ctx.info(message, linter=cls.name(), node=requirements_nodes)


class CondaRequirementInexact(Linter):
    @classmethod
    def lint(cls, tool_source: "ToolSource", lint_ctx: "LintContext"):
        for conda_target in _requirements_conda_targets(tool_source):
            (best_hit, exact) = best_practice_search(conda_target)
            conda_target_str = conda_target.package
            if conda_target.version:
                conda_target_str += "@%s" % (conda_target.version)
            if best_hit and not exact:
                message = f"Requirement [{conda_target_str}] doesn't exactly match available version [{best_hit['version']}] in best practice Conda channel [{best_hit.get('channel')}]."
                requirements_node = xml_node_from_toolsource(tool_source, "requirements")
                lint_ctx.warn(message, linter=cls.name(), node=requirements_node)


class CondaRequirementMissing(Linter):
    @classmethod
    def lint(cls, tool_source: "ToolSource", lint_ctx: "LintContext"):
        for conda_target in _requirements_conda_targets(tool_source):
            (best_hit, exact) = best_practice_search(conda_target)
            conda_target_str = conda_target.package
            if conda_target.version:
                conda_target_str += "@%s" % (conda_target.version)
            if best_hit and not exact:
                message = f"Requirement [{conda_target_str}] doesn't match any recipe in a best practice conda channel ['{BEST_PRACTICE_CHANNELS}']."
                requirements_node = xml_node_from_toolsource(tool_source, "requirements")
                lint_ctx.warn(message, linter=cls.name(), node=requirements_node)


def _requirements_conda_targets(tool_source: "ToolSource") -> Generator[CondaTarget]:
    requirements, *_ = tool_source.parse_requirements_and_containers()
    for requirement in requirements:
        conda_target = requirement_to_conda_targets(requirement)
        if conda_target:
            yield conda_target
