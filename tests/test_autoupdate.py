from planemo.autoupdate import get_newest_tool_id


def test_get_newest_tool_id():
    tool_ids = [
        "toolshed.g2.bx.psu.edu/repos/iuc/rgrnastar/rna_star/2.7.8a+galaxy0",
        "toolshed.g2.bx.psu.edu/repos/iuc/rgrnastar/rna_star/2.7.8a",
    ]
    assert get_newest_tool_id(tool_ids) == tool_ids[0]
