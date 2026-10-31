"""Lint Galaxy tool help text URLs."""

from typing import (
    Any,
    List,
    NamedTuple,
    Optional,
    TYPE_CHECKING,
)

from galaxy.tool_util.lint import Linter

from planemo.lint import (
    BROWSER_USER_AGENT,
)
from planemo.lint import validate_url as _validate_url
from planemo.linters.util import cached_lint_result
from planemo.shed import _find_urls_in_text

if TYPE_CHECKING:
    from galaxy.tool_util.lint import LintContext
    from galaxy.tool_util.parser.interface import ToolSource


class URLCheck(NamedTuple):
    url: str
    help_node: Any
    error: Optional[str]


def _url_checks(tool_source: "ToolSource") -> List[URLCheck]:
    def calculate() -> List[URLCheck]:
        tool_xml = getattr(tool_source, "xml_tree", None)
        if not tool_xml:
            return []

        checks = []
        for help_node in tool_xml.findall("help"):
            for url_match in _find_urls_in_text(help_node.text or ""):
                url = url_match[0]
                checks.append(URLCheck(url, help_node, _validate_url(url, BROWSER_USER_AGENT)))
        return checks

    return cached_lint_result(tool_source, "urls", calculate)


class URLInaccessible(Linter):
    @classmethod
    def lint(cls, tool_source: "ToolSource", lint_ctx: "LintContext"):
        for check in _url_checks(tool_source):
            if check.error:
                lint_ctx.error(check.error, linter=cls.name(), node=check.help_node)


class URLValid(Linter):
    @classmethod
    def lint(cls, tool_source: "ToolSource", lint_ctx: "LintContext"):
        for check in _url_checks(tool_source):
            if not check.error:
                lint_ctx.info(f"URL OK {check.url}", linter=cls.name(), node=check.help_node)
