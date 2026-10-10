from typing import (
    Any,
    Callable,
    Dict,
    Optional,
    TYPE_CHECKING,
    TypeVar,
)
from weakref import WeakKeyDictionary

if TYPE_CHECKING:
    from galaxy.tool_util.parser.interface import ToolSource
    from galaxy.util import Element

LintResult = TypeVar("LintResult")
_lint_result_cache: "WeakKeyDictionary[ToolSource, Dict[str, Any]]" = WeakKeyDictionary()


def cached_lint_result(tool_source: "ToolSource", key: str, calculate: Callable[[], LintResult]) -> LintResult:
    """Cache shared analysis while a tool source remains alive."""
    tool_cache = _lint_result_cache.setdefault(tool_source, {})
    if key not in tool_cache:
        tool_cache[key] = calculate()
    return tool_cache[key]


def xml_node_from_toolsource(tool_source: "ToolSource", tag: str) -> Optional["Element"]:
    node = None
    xml_tree = getattr(tool_source, "xml_tree", None)
    if xml_tree:
        node = xml_tree.find(tag)
    return node
