"""Tests for the changed-path selection in :mod:`planemo.ci`."""

import os

import pytest

from planemo.ci import (
    changed_repos,
    changed_repos_extended,
    changed_tools_extended,
)

TOOL_XML = """<tool id="{id}" name="{id}" version="1.0.0">
    <command>echo hello</command>
    <inputs />
    <outputs />
</tool>
"""

BROKEN_TOOL_XML = """<tool id="broken" name="broken" version="1.0.0">
    <command>echo hello</command>
    <inputs>
</tool>
"""


def _write(path, contents=""):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(contents)


def _shed_repo(root, rel, tool_id, extra_files=()):
    """Write a shed repository with a single tool at ``rel`` under ``root``."""
    _write(os.path.join(root, rel, ".shed.yml"), "name: %s\n" % tool_id)
    _write(os.path.join(root, rel, "%s.xml" % tool_id), TOOL_XML.format(id=tool_id))
    for extra in extra_files:
        _write(os.path.join(root, rel, extra), "data\n")


@pytest.fixture
def repo(tmp_path, monkeypatch):
    """A small tools-iuc shaped tree, with cwd at its root."""
    root = str(tmp_path)
    # A repository with no subdirectories at all.
    _shed_repo(root, "tools/flat", "flat")
    # A repository owning a test-data directory and a helper script.
    _shed_repo(root, "tools/withdata", "withdata", extra_files=["test-data/in.txt", "scripts/helper.py"])
    # A tool collection: shared macros in the parent, .shed.yml per sub-tool.
    _write(os.path.join(root, "tool_collections/coll/macros.xml"), "<macros />\n")
    _shed_repo(root, "tool_collections/coll/one", "one")
    _shed_repo(root, "tool_collections/coll/two", "two")
    monkeypatch.chdir(root)
    return root


def test_changed_repos_default_finds_owning_repo(repo):
    """Default behavior walks up to the owning .shed.yml -- master's semantics."""
    assert changed_repos(["tools/flat/flat.xml"]) == {"tools/flat"}
    assert changed_repos(["tools/withdata/test-data/in.txt"]) == {"tools/withdata"}


def test_changed_repos_default_misses_sibling_repos(repo):
    """Default behavior does not descend -- this is the gap --extended_git_diff closes."""
    assert changed_repos(["tool_collections/coll/macros.xml"]) == set()


def test_changed_repos_extended_returns_normalized_paths(repo):
    """Regression: a repo with no subdirectory must not come back slash-suffixed.

    ``filter_paths`` intersects against ``os.path.relpath`` output, so a trailing
    separator silently drops the repository from the selection.
    """
    assert changed_repos_extended(["tools/flat/flat.xml"]) == {"tools/flat"}


def test_changed_repos_extended_descends_into_subdirectories(repo):
    """A macros change in the collection parent selects each sub-repository."""
    expected = {"tool_collections/coll/one", "tool_collections/coll/two"}
    assert changed_repos_extended(["tool_collections/coll/macros.xml"]) == expected


def test_changed_tools_extended_resolves_non_tool_file_to_owning_tool(repo):
    """Issue #1129: a test-data or script change should select the owning tool."""
    cwd = os.getcwd()
    assert changed_tools_extended(None, ["tools/withdata/test-data/in.txt"], cwd) == {"tools/withdata/withdata.xml"}
    assert changed_tools_extended(None, ["tools/withdata/scripts/helper.py"], cwd) == {"tools/withdata/withdata.xml"}


def test_changed_tools_extended_survives_deleted_directory(repo):
    """``git diff`` lists deleted files, so the directory may be gone."""
    assert changed_tools_extended(None, ["tools/removed/removed.xml"], os.getcwd()) == set()


def test_changed_tools_extended_does_not_escalate_on_unparsable_tool(repo, tmp_path):
    """A tool broken by the commit under test must not select its siblings.

    An unparsable XML makes the directory look empty to the walk, which would
    otherwise escalate to the parent and select every unrelated tool under it.
    """
    _write(os.path.join(str(tmp_path), "tools/broken/broken.xml"), BROKEN_TOOL_XML)
    selected = changed_tools_extended(None, ["tools/broken/broken.xml"], os.getcwd())
    assert "tools/flat/flat.xml" not in selected
    assert "tools/withdata/withdata.xml" not in selected


def test_changed_tools_extended_deduplicates_directories(repo):
    """Several files in one directory resolve to the same tool exactly once."""
    changed = [
        "tools/withdata/test-data/in.txt",
        "tools/withdata/scripts/helper.py",
        "tools/withdata/.shed.yml",
    ]
    assert changed_tools_extended(None, changed, os.getcwd()) == {"tools/withdata/withdata.xml"}
