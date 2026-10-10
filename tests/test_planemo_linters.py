from pathlib import Path
from types import SimpleNamespace
from unittest import mock

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
from planemo.shed_lint import (
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
