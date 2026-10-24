import pytest

from planemo.tools import (
    _parse_tool_version,
    parse_tool_version,
)

ORDERED = ["1.0", "1.0+galaxy0", "1.0+galaxy1", "1.0+galaxy9", "1.0+galaxy10", "1.0+galaxy11", "1.1"]


@pytest.mark.parametrize("parse", [parse_tool_version, _parse_tool_version])
def test_galaxy_suffix_orders_numerically(parse):
    assert sorted(reversed(ORDERED), key=parse) == ORDERED
    assert parse("1.0+galaxy10") > parse("1.0+galaxy9")
