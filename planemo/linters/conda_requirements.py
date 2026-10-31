"""Ensure requirements are matched in best-practice Conda channels."""

from typing import (
    Any,
    List,
    NamedTuple,
    Optional,
    TYPE_CHECKING,
)

from galaxy.tool_util.lint import Linter

from planemo.conda import (
    BEST_PRACTICE_CHANNELS,
    best_practice_search,
    tool_source_conda_targets,
)
from planemo.linters.util import (
    cached_lint_result,
    xml_node_from_toolsource,
)

if TYPE_CHECKING:
    from galaxy.tool_util.deps.conda_util import CondaTarget
    from galaxy.tool_util.lint import LintContext
    from galaxy.tool_util.parser.interface import ToolSource

MESSAGE_WARN_NO_REQUIREMENTS = "No valid package requirement tags found to check against Conda."

lint_tool_types = ["*"]


class CondaCheck(NamedTuple):
    target: "CondaTarget"
    best_hit: Optional[Any]
    exact: Optional[bool]


def _conda_checks(tool_source: "ToolSource") -> List[CondaCheck]:
    def calculate() -> List[CondaCheck]:
        checks = []
        for conda_target in tool_source_conda_targets(tool_source):
            best_hit, exact = best_practice_search(conda_target)
            checks.append(CondaCheck(conda_target, best_hit, exact))
        return checks

    return cached_lint_result(tool_source, "conda_requirements", calculate)


def _target_string(check: CondaCheck) -> str:
    target = check.target.package
    if check.target.version:
        target += f"@{check.target.version}"
    return target


class CondaRequirementInexact(Linter):
    @classmethod
    def lint(cls, tool_source: "ToolSource", lint_ctx: "LintContext"):
        requirements_node = xml_node_from_toolsource(tool_source, "requirements")
        for check in _conda_checks(tool_source):
            if check.best_hit and not check.exact:
                message = (
                    f"Requirement [{_target_string(check)}] doesn't exactly match available version "
                    f"[{check.best_hit['version']}] in best practice Conda channel "
                    f"[{check.best_hit.get('channel')}]."
                )
                lint_ctx.warn(message, linter=cls.name(), node=requirements_node)


class CondaRequirementMissing(Linter):
    @classmethod
    def lint(cls, tool_source: "ToolSource", lint_ctx: "LintContext"):
        requirements_node = xml_node_from_toolsource(tool_source, "requirements")
        for check in _conda_checks(tool_source):
            if not check.best_hit:
                message = (
                    f"Requirement [{_target_string(check)}] doesn't match any recipe in a best practice "
                    f"Conda channel [{BEST_PRACTICE_CHANNELS}]."
                )
                lint_ctx.warn(message, linter=cls.name(), node=requirements_node)


class CondaRequirementValid(Linter):
    @classmethod
    def lint(cls, tool_source: "ToolSource", lint_ctx: "LintContext"):
        requirements_node = xml_node_from_toolsource(tool_source, "requirements")
        for check in _conda_checks(tool_source):
            if check.best_hit and check.exact:
                message = (
                    f"Requirement [{_target_string(check)}] matches target in best practice "
                    f"Conda channel [{check.best_hit.get('channel')}]."
                )
                lint_ctx.info(message, linter=cls.name(), node=requirements_node)


class CondaRequirementsMissing(Linter):
    @classmethod
    def lint(cls, tool_source: "ToolSource", lint_ctx: "LintContext"):
        if not _conda_checks(tool_source):
            requirements_node = xml_node_from_toolsource(tool_source, "requirements")
            lint_ctx.warn(MESSAGE_WARN_NO_REQUIREMENTS, linter=cls.name(), node=requirements_node)
