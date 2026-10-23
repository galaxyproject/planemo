"""Unit tests for the ShedVersion check's comparison against the published tool."""

import responses

from planemo.shed_lint import (
    _published_tool_source,
    _same_tool,
)

TOOL = '<tool id="cat" name="Cat" version="1.0"><command>cat $input &gt; $output</command></tool>'


class _Tsi:
    base_url = "https://toolshed.example.org/"


def test_same_tool_ignores_formatting_whitespace():
    reformatted = '<tool id="cat" name="Cat" version="1.0">\n    <command>cat $input &gt; $output</command>\n</tool>\n'
    assert _same_tool(TOOL, reformatted)


def test_same_tool_detects_a_change():
    assert not _same_tool(TOOL, TOOL.replace("cat $input", "cat -n $input"))


def test_same_tool_treats_unparseable_source_as_changed():
    assert not _same_tool(TOOL, "<html>not a tool")
    assert not _same_tool("<tool", TOOL)


@responses.activate
def test_published_tool_source_fetches_by_trs_id_and_version():
    url = "https://toolshed.example.org/api/tools/owner~repo~cat/versions/1.0%2Bgalaxy1/tool_source"
    responses.add(responses.GET, url, body=TOOL, status=200)
    assert _published_tool_source(_Tsi(), "owner", "repo", "cat", "1.0+galaxy1") == TOOL


@responses.activate
def test_published_tool_source_is_none_when_unavailable():
    url = "https://toolshed.example.org/api/tools/owner~repo~cat/versions/1.0/tool_source"
    responses.add(responses.GET, url, status=404)
    assert _published_tool_source(_Tsi(), "owner", "repo", "cat", "1.0") is None
