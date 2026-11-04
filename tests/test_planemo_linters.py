from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import pytest
import responses
from galaxy.tool_util.lint import (
    lint_tool_source_with_modules,
    LintContext,
)
from galaxy.tool_util.parser import get_tool_source

from planemo.lint import build_lint_args
from planemo.linters import (
    biocontainer_registered,
    conda_requirements,
    doi,
    urls,
)
from planemo.shed.interface import (
    shed_tool_source,
    tool_shed_instance,
)
from planemo.shed_lint import (
    lint_shed_version,
    lint_tool_dependencies_urls,
    REPOSITORY_LINTER_NAMES,
)
from .test_utils import (
    PROJECT_TEMPLATES_DIR,
    TEST_REPOS_DIR,
    TEST_TOOLS_DIR,
)


def _tool_source(directory, name):
    return get_tool_source(str(Path(directory) / name))


def _messages(lint_ctx):
    return [(message.level, message.message, message.linter) for message in lint_ctx.message_list]


def test_doi_linter_checks_each_citation_once():
    tool_source = _tool_source(TEST_TOOLS_DIR, "invalid_doi.xml")
    responses = [mock.Mock(status_code=200), mock.Mock(status_code=404)]
    lint_ctx = LintContext("all")

    with mock.patch("planemo.linters.doi.requests.get", side_effect=responses) as get:
        lint_tool_source_with_modules(lint_ctx, tool_source, [doi])

    assert get.call_count == 2
    assert all(call.kwargs["timeout"] == 5 for call in get.call_args_list)
    assert _messages(lint_ctx) == [
        ("error", "doi:10.1101/666666 is not a valid DOI", "DoiInvalid"),
        (
            "error",
            "doi:10.1101/014043 is valid, but Galaxy expects DOI without 'doi:' prefix",
            "DoiPrefix",
        ),
    ]


def test_doi_linter_can_be_skipped_by_module_name():
    tool_source = _tool_source(TEST_TOOLS_DIR, "invalid_doi.xml")
    lint_ctx = LintContext("all", skip_types=["doi"])

    with mock.patch("planemo.linters.doi.requests.get") as get:
        lint_tool_source_with_modules(lint_ctx, tool_source, [doi])

    get.assert_not_called()
    assert lint_ctx.message_list == []


def test_doi_linter_reports_request_failures_without_crashing():
    tool_source = _tool_source(TEST_TOOLS_DIR, "invalid_doi.xml")
    lint_ctx = LintContext("all")

    with mock.patch(
        "planemo.linters.doi.requests.get",
        side_effect=doi.requests.ConnectionError("offline"),
    ) as get:
        lint_tool_source_with_modules(lint_ctx, tool_source, [doi])

    assert get.call_count == 2
    assert _messages(lint_ctx) == [
        (
            "warning",
            "Error 'offline' accessing https://doi.org/10.1101/014043",
            "DoiUnexpectedResponse",
        ),
        (
            "warning",
            "Error 'offline' accessing https://doi.org/10.1101/666666",
            "DoiUnexpectedResponse",
        ),
    ]


def test_url_linter_checks_each_help_url_once():
    tool_source = _tool_source(TEST_TOOLS_DIR, "url.xml")
    response = mock.Mock(status_code=200)
    response.iter_content.side_effect = lambda _: iter([b"content"])
    lint_ctx = LintContext("all")

    with mock.patch("planemo.lint.requests.get", return_value=response) as get:
        lint_tool_source_with_modules(lint_ctx, tool_source, [urls])

    assert get.call_count == 2
    assert all(call.kwargs["timeout"] == 5 for call in get.call_args_list)
    assert _messages(lint_ctx) == [
        ("info", "URL OK http://galaxyproject.org/", "URLValid"),
        ("info", "URL OK https://galaxyproject.org/", "URLValid"),
    ]


def test_tool_dependency_url_linter_reuses_hardened_url_check():
    realized_repository = SimpleNamespace(real_path=str(Path(TEST_REPOS_DIR) / "package_1"))
    response = mock.Mock(status_code=200)
    response.iter_content.side_effect = lambda _: iter([b"content"])
    lint_ctx = LintContext("all")

    with mock.patch("planemo.lint.requests.get", return_value=response) as get:
        lint_tool_dependencies_urls(realized_repository, lint_ctx)

    assert get.call_count == 5
    assert all(call.kwargs["timeout"] == 5 for call in get.call_args_list)
    assert all(message.level == "info" for message in lint_ctx.message_list)


def test_repository_linter_names_are_valid_skip_targets():
    ctx = mock.Mock(global_config={})

    with mock.patch("planemo.lint.error") as error:
        lint_args = build_lint_args(
            ctx,
            skip=["shed_remote_repository_url"],
            extra_linter_names=REPOSITORY_LINTER_NAMES,
        )

    error.assert_not_called()
    assert lint_args["skip_types"] == ["shed_remote_repository_url"]


def test_conda_linter_checks_each_requirement_once():
    tool_source = _tool_source(TEST_TOOLS_DIR, "bwa_invalid_version.xml")
    lint_ctx = LintContext("all")

    with mock.patch("planemo.linters.conda_requirements.best_practice_search", return_value=(None, None)) as search:
        lint_tool_source_with_modules(lint_ctx, tool_source, [conda_requirements])

    search.assert_called_once()
    assert _messages(lint_ctx) == [
        (
            "warning",
            "Requirement [bwa@0.4.12] doesn't match any recipe in a best practice "
            "Conda channel [['conda-forge', 'bioconda']].",
            "CondaRequirementMissing",
        )
    ]


def test_conda_linter_distinguishes_inexact_from_missing_requirements():
    tool_source = _tool_source(TEST_TOOLS_DIR, "bwa_invalid_version.xml")
    lint_ctx = LintContext("all")
    best_hit = {"channel": "bioconda", "version": "0.7.10"}

    with mock.patch(
        "planemo.linters.conda_requirements.best_practice_search",
        return_value=(best_hit, False),
    ) as search:
        lint_tool_source_with_modules(lint_ctx, tool_source, [conda_requirements])

    search.assert_called_once()
    assert _messages(lint_ctx) == [
        (
            "warning",
            "Requirement [bwa@0.4.12] doesn't exactly match available version "
            "[0.7.10] in best practice Conda channel [bioconda].",
            "CondaRequirementInexact",
        )
    ]


def test_conda_linter_reports_missing_requirements_without_searching():
    tool_source = _tool_source(TEST_TOOLS_DIR, "bwa_without_requirements.xml")
    lint_ctx = LintContext("all")

    with mock.patch("planemo.linters.conda_requirements.best_practice_search") as search:
        lint_tool_source_with_modules(lint_ctx, tool_source, [conda_requirements])

    search.assert_not_called()
    assert _messages(lint_ctx) == [
        (
            "warning",
            "No valid package requirement tags found to check against Conda.",
            "CondaRequirementsMissing",
        )
    ]


def test_biocontainer_linter_resolves_requirement_set_once():
    tool_source = _tool_source(Path(PROJECT_TEMPLATES_DIR) / "seqtk_complete", "seqtk_seq.xml")
    lint_ctx = LintContext("all")

    with mock.patch(
        "planemo.linters.biocontainer_registered.mulled_container_name",
        return_value="quay.io/biocontainers/seqtk:1.2",
    ) as resolve:
        lint_tool_source_with_modules(lint_ctx, tool_source, [biocontainer_registered])

    resolve.assert_called_once()
    assert _messages(lint_ctx) == [
        (
            "info",
            "BioContainer best-practice container found [quay.io/biocontainers/seqtk:1.2].",
            "BiocontainerValid",
        )
    ]


def test_biocontainer_linter_reports_missing_requirements_without_resolving():
    tool_source = _tool_source(TEST_TOOLS_DIR, "bwa_without_requirements.xml")
    lint_ctx = LintContext("all")

    with mock.patch("planemo.linters.biocontainer_registered.mulled_container_name") as resolve:
        lint_tool_source_with_modules(lint_ctx, tool_source, [biocontainer_registered])

    resolve.assert_not_called()
    assert _messages(lint_ctx) == [
        (
            "warning",
            "No valid package requirement tags found to infer BioContainer from.",
            "BiocontainerRequirementsMissing",
        )
    ]


def test_biocontainer_linter_reports_missing_container_once():
    tool_source = _tool_source(Path(PROJECT_TEMPLATES_DIR) / "seqtk_complete", "seqtk_seq.xml")
    lint_ctx = LintContext("all")

    with mock.patch(
        "planemo.linters.biocontainer_registered.mulled_container_name",
        return_value=None,
    ) as resolve:
        lint_tool_source_with_modules(lint_ctx, tool_source, [biocontainer_registered])

    resolve.assert_called_once()
    assert _messages(lint_ctx) == [
        (
            "warning",
            "Failed to find a BioContainer registered for these requirements.",
            "BiocontainerMissing",
        )
    ]


SHED_VERSION_TOOL = """<tool id="shed_version_tool" name="Shed Version Tool" version="1.0+galaxy10">
    <command>echo hi</command>
    <inputs />
    <outputs />
</tool>
"""


@dataclass
class _ShedRepository:
    path: str
    owner: str = "iuc"
    name: str = "shed_version"


SHED_VERSION_OTHER_TOOL = """<tool id="shed_version_other_tool" name="Shed Version Other Tool" version="2.0">
    <command>echo bye</command>
    <inputs />
    <outputs />
</tool>
"""


class _FakeShedRepositories:
    def __init__(self, valid_tools):
        self.valid_tools = valid_tools

    def get_repositories(self, name, owner):
        return [{"name": name, "owner": owner}]

    def get_ordered_installable_revisions(self, name, owner):
        return ["abc123"]

    def get_repository_revision_install_info(self, name, owner, revision):
        return {}, {"valid_tools": [dict(tool) for tool in self.valid_tools]}, {}


class _FakeToolShed:
    def __init__(self, valid_tools=None):
        if valid_tools is None:
            valid_tools = [{"id": "shed_version_tool", "version": "1.0+galaxy10"}]
        self.repositories = _FakeShedRepositories(valid_tools)


def _lint_shed_version(tmp_path, tool_xmls=None, ts_tool_source=None, valid_tools=None):
    tool_xmls = tool_xmls if tool_xmls is not None else {"shed_version_tool.xml": SHED_VERSION_TOOL}
    for name, tool_xml in tool_xmls.items():
        (tmp_path / name).write_text(tool_xml)
    lint_ctx = LintContext("all")

    with (
        mock.patch("planemo.shed_lint.tool_shed_instance", return_value=_FakeToolShed(valid_tools)),
        mock.patch("planemo.shed_lint.shed_tool_source", return_value=ts_tool_source) as shed_source,
    ):
        lint_ctx.lint("shed_version", lint_shed_version, _ShedRepository(str(tmp_path)))

    return lint_ctx, shed_source


@pytest.mark.parametrize(
    "ts_version, expected",
    [
        ("1.0+galaxy9", []),
        (
            "1.0+galaxy10",
            [
                (
                    "error",
                    "shed_version_tool: version 1.0+galaxy10 is less or equal than version of the latest installable revision 1.0+galaxy10",
                    "ShedVersion",
                )
            ],
        ),
    ],
)
def test_shed_version_compares_galaxy_suffix_numerically(tmp_path, ts_version, expected):
    # the published source differs from the local one, so a bump is required
    lint_ctx, shed_source = _lint_shed_version(
        tmp_path,
        ts_tool_source='<tool id="shed_version_tool" name="Shed Version Tool" version="1.0+galaxy10"/>',
        valid_tools=[{"id": "shed_version_tool", "version": ts_version}],
    )

    assert _messages(lint_ctx) == expected
    if expected == []:
        # a bumped version needs no comparison with the published source
        shed_source.assert_not_called()


def test_shed_version_unchanged_tool_source_needs_no_bump(tmp_path):
    # e.g. another tool of a multi tool repository was bumped: this tool was not
    # touched, so it must be neither bumped nor reported.
    (tmp_path / "shed_version_tool.xml").write_text(SHED_VERSION_TOOL)
    published_source = get_tool_source(str(tmp_path / "shed_version_tool.xml")).to_string()
    lint_ctx, shed_source = _lint_shed_version(tmp_path, ts_tool_source=published_source)

    assert _messages(lint_ctx) == [
        (
            "info",
            "shed_version_tool: tool source unchanged since installable revision 1.0+galaxy10, no version bump needed.",
            "ShedVersion",
        )
    ]
    shed_source.assert_called_once()


def test_shed_version_changed_tool_source_needs_bump(tmp_path):
    lint_ctx, _ = _lint_shed_version(
        tmp_path, ts_tool_source='<tool id="shed_version_tool" name="Shed Version Tool" version="1.0+galaxy10"/>'
    )

    assert _messages(lint_ctx) == [
        (
            "error",
            "shed_version_tool: version 1.0+galaxy10 is less or equal than version of the latest installable revision 1.0+galaxy10",
            "ShedVersion",
        )
    ]


def test_shed_version_unfetchable_tool_source_is_an_error(tmp_path):
    lint_ctx, _ = _lint_shed_version(tmp_path, ts_tool_source=None)

    assert _messages(lint_ctx) == [
        (
            "error",
            "shed_version_tool: could not fetch the tool source published with version 1.0+galaxy10 "
            "from the Tool Shed, cannot check whether a version bump is needed.",
            "ShedVersion",
        )
    ]


def test_shed_version_untouched_tools_of_multi_tool_repository_are_clean(tmp_path):
    (tmp_path / "shed_version_tool.xml").write_text(SHED_VERSION_TOOL)
    published_source = get_tool_source(str(tmp_path / "shed_version_tool.xml")).to_string()

    lint_ctx, _ = _lint_shed_version(
        tmp_path,
        tool_xmls={
            "shed_version_tool.xml": SHED_VERSION_TOOL,
            # this tool is bumped and must not be reported
            "shed_version_other_tool.xml": SHED_VERSION_OTHER_TOOL,
        },
        ts_tool_source=published_source,
        valid_tools=[
            {"id": "shed_version_tool", "version": "1.0+galaxy10"},
            {"id": "shed_version_other_tool", "version": "1.0"},
        ],
    )

    assert _messages(lint_ctx) == [
        (
            "info",
            "shed_version_tool: tool source unchanged since installable revision 1.0+galaxy10, no version bump needed.",
            "ShedVersion",
        )
    ]


def test_shed_version_new_tool_in_existing_repository_is_clean(tmp_path):
    # the previous implementation asserted here, because the tool is not part of
    # the latest installable revision of the repository yet
    lint_ctx, shed_source = _lint_shed_version(tmp_path, valid_tools=[{"id": "some_other_tool", "version": "1.0"}])

    assert _messages(lint_ctx) == [
        ("info", "shed_version_tool: not yet published in shed_version, no version check needed.", "ShedVersion")
    ]
    shed_source.assert_not_called()


def test_shed_version_downgrade_is_an_error_without_source_comparison(tmp_path):
    lint_ctx, shed_source = _lint_shed_version(
        tmp_path, tool_xmls={"shed_version_tool.xml": SHED_VERSION_TOOL.replace("1.0+galaxy10", "1.0+galaxy9")}
    )

    assert _messages(lint_ctx) == [
        (
            "error",
            "shed_version_tool: version 1.0+galaxy9 is less or equal than version of the latest installable revision 1.0+galaxy10",
            "ShedVersion",
        )
    ]
    shed_source.assert_not_called()


def test_shed_tool_source_requests_published_source():
    tsi = tool_shed_instance("https://toolshed.g2.bx.psu.edu")
    url = (
        "https://toolshed.g2.bx.psu.edu/api/tools/iuc~shed_version~shed_version_tool"
        "/versions/1.0%2Bgalaxy10/tool_source"
    )

    with responses.RequestsMock() as rsps:
        rsps.add(responses.GET, url, body="<tool/>", status=200, headers={"language": "xml"})
        assert shed_tool_source(tsi, "iuc", "shed_version", "shed_version_tool", "1.0+galaxy10") == "<tool/>"

    # e.g. a Tool Shed that does not provide the endpoint (< 26.1) or unknown tools
    with responses.RequestsMock() as rsps:
        rsps.add(responses.GET, url, status=404)
        assert shed_tool_source(tsi, "iuc", "shed_version", "shed_version_tool", "1.0+galaxy10") is None
