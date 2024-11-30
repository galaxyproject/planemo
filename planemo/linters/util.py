from typing import (
    Optional,
    TYPE_CHECKING,
)

if TYPE_CHECKING:
    from galaxy.tool_util.parser.interface import ToolSource
    from galaxy.util import Element


def xml_node_from_toolsource(tool_source: "ToolSource", tag: "str") -> Optional["Element"]:
    node = None
    xml_tree = getattr(tool_source, "xml_tree", None)
    if xml_tree:
        node = xml_tree.find(tag)
    return node
