from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import pytest
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
    tool_data_tables,
    urls,
)
from planemo.shed_lint import (
    lint_shed_version,
    lint_tool_dependencies_urls,
    REPOSITORY_LINTER_NAMES,
)
from planemo.tool_lint import _lint_extra_modules
from .test_utils import (
    PROJECT_TEMPLATES_DIR,
    TEST_DATA_DIR,
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


class _FakeShedRepositories:
    def __init__(self, ts_version):
        self.ts_version = ts_version

    def get_repositories(self, name, owner):
        return [{"name": name, "owner": owner}]

    def get_ordered_installable_revisions(self, name, owner):
        return ["abc123"]

    def get_repository_revision_install_info(self, name, owner, revision):
        return {}, {"valid_tools": [{"id": "shed_version_tool", "version": self.ts_version}]}, {}


class _FakeToolShed:
    def __init__(self, ts_version):
        self.repositories = _FakeShedRepositories(ts_version)


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
    (tmp_path / "shed_version_tool.xml").write_text(SHED_VERSION_TOOL)
    lint_ctx = LintContext("all")

    with mock.patch("planemo.shed_lint.tool_shed_instance", return_value=_FakeToolShed(ts_version)):
        lint_ctx.lint("shed_version", lint_shed_version, _ShedRepository(str(tmp_path)))

    assert _messages(lint_ctx) == expected


DATA_TABLE_FIXTURES_DIR = Path(TEST_DATA_DIR) / "lint_data_tables"


def _lint_data_table_tool(case, name):
    tool_source = _tool_source(DATA_TABLE_FIXTURES_DIR / case, name)
    lint_ctx = LintContext("all")
    lint_tool_source_with_modules(lint_ctx, tool_source, [tool_data_tables])
    return lint_ctx


def _lint_levels_linters_messages(lint_ctx):
    return [(level, linter, message) for level, message, linter in _messages(lint_ctx)]


def test_data_table_linter_accepts_repository_with_sample_files():
    lint_ctx = _lint_data_table_tool("ok_data_table", "gene_ontology.xml")

    assert lint_ctx.message_list == []


def test_data_table_linter_requires_data_table_config_sample():
    lint_ctx = _lint_data_table_tool("missing_conf", "gene_ontology.xml")

    [(level, linter, message)] = _lint_levels_linters_messages(lint_ctx)
    assert (level, linter) == ("error", "ToolDataTableConfSampleMissing")
    assert "tool_data_table_conf.xml.sample" in message
    assert "[gene_ontology]" in message


def test_data_table_linter_requires_index_file_sample():
    lint_ctx = _lint_data_table_tool("missing_loc_sample", "gene_ontology.xml")

    [(level, linter, message)] = _lint_levels_linters_messages(lint_ctx)
    assert (level, linter) == ("error", "ToolDataTableIndexSampleMissing")
    assert "gene_ontology.loc.sample" in message


def test_data_table_linter_reports_unparsable_config_sample():
    lint_ctx = _lint_data_table_tool("invalid_conf", "gene_ontology.xml")

    [(level, linter, message)] = _lint_levels_linters_messages(lint_ctx)
    assert (level, linter) == ("error", "ToolDataTableConfSampleInvalid")
    assert "tool_data_table_conf.xml.sample" in message


def test_data_table_linter_warns_about_config_naming_sample_file():
    lint_ctx = _lint_data_table_tool("conf_refers_to_sample", "gene_ontology.xml")

    [(level, linter, message)] = _lint_levels_linters_messages(lint_ctx)
    assert (level, linter) == ("warning", "ToolDataTableConfRefersToSample")
    assert "names the sample file gene_ontology.loc.sample instead of the data file" in message
    assert "installing this repository creates gene_ontology.loc" in message


def test_data_table_linter_warns_about_from_file_without_sample():
    lint_ctx = _lint_data_table_tool("from_file_missing_sample", "from_file.xml")

    [(level, linter, message)] = _lint_levels_linters_messages(lint_ctx)
    assert (level, linter) == ("warning", "ToolDataTableFromFileSampleMissing")
    assert "builds.loc.sample" in message


def test_data_table_linter_ignores_tables_of_galaxy_or_other_repositories():
    # Tables declared by Galaxy itself or by another repository need no entry
    # in the data table configuration of this repository and the Tool Shed does
    # not require sample files for them either.
    lint_ctx = _lint_data_table_tool("table_from_other_repo", "all_fasta.xml")

    assert lint_ctx.message_list == []


def test_data_table_linter_can_be_skipped_by_module_name():
    tool_source = _tool_source(DATA_TABLE_FIXTURES_DIR / "missing_loc_sample", "gene_ontology.xml")
    lint_ctx = LintContext("all", skip_types=["tool_data_tables"])

    lint_tool_source_with_modules(lint_ctx, tool_source, [tool_data_tables])

    assert lint_ctx.message_list == []


def test_data_table_linter_skip_names_are_valid():
    ctx = mock.Mock(global_config={})

    with mock.patch("planemo.lint.error") as error:
        lint_args = build_lint_args(
            ctx,
            skip=["tool_data_tables", "ToolDataTableIndexSampleMissing", "ToolDataTableConfRefersToSample"],
        )

    error.assert_not_called()
    assert "ToolDataTableIndexSampleMissing" in lint_args["skip_types"]


def test_data_table_linter_is_enabled_by_default():
    assert tool_data_tables in _lint_extra_modules()
