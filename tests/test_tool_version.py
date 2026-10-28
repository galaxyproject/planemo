import pytest
from galaxy.tool_util import version as galaxy_tool_util_version

from planemo.tools import _parse_tool_version

ORDERED = ["1.0", "1.0+galaxy0", "1.0+galaxy1", "1.0+galaxy9", "1.0+galaxy10", "1.0+galaxy11", "1.1"]


def test_fallback_orders_galaxy_suffix_numerically():
    assert sorted(reversed(ORDERED), key=_parse_tool_version) == ORDERED


def test_fallback_matches_upstream_parse_tool_version():
    upstream = getattr(galaxy_tool_util_version, "parse_tool_version", None)
    if upstream is None:
        pytest.skip("galaxy-tool-util lacks parse_tool_version")
    assert sorted(reversed(ORDERED), key=upstream) == ORDERED
