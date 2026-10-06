"""Tests for runnable type detection."""

# Planemo's normal startup loads the Galaxy package before importing runnable.
from planemo.galaxy import galaxy_config  # noqa: F401
from planemo.runnable import (
    for_path,
    RunnableType,
)


def test_yaml_galaxy_tool_is_a_runnable(tmp_path):
    tool_path = tmp_path / "minimal.yml"
    tool_path.write_text("class: GalaxyTool\nid: minimal\nname: Minimal\nversion: '1.0'\n")

    assert for_path(str(tool_path)).type == RunnableType.galaxy_tool
