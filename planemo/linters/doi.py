"""Lint DOI citations in Galaxy tools."""

from typing import (
    Any,
    List,
    NamedTuple,
    Optional,
    TYPE_CHECKING,
)

import requests
from galaxy.tool_util.lint import Linter

from planemo.lint import REQUEST_TIMEOUT
from planemo.linters.util import cached_lint_result

if TYPE_CHECKING:
    from galaxy.tool_util.lint import LintContext
    from galaxy.tool_util.parser.interface import ToolSource


class DoiCheck(NamedTuple):
    citation: Any
    publication_id: str
    doiless_publication_id: str
    status_code: Optional[int]
    error: Optional[str]


def _doi_checks(tool_source: "ToolSource") -> List[DoiCheck]:
    def calculate() -> List[DoiCheck]:
        tool_xml = getattr(tool_source, "xml_tree", None)
        if not tool_xml:
            return []
        citations = tool_xml.find("citations")
        if citations is None:
            return []

        checks = []
        for citation in citations:
            if citation.tag != "citation" or citation.attrib.get("type") != "doi":
                continue
            publication_id = (citation.text or "").strip()
            doiless_publication_id = publication_id.split("doi:", 1)[-1]
            if not doiless_publication_id:
                checks.append(DoiCheck(citation, publication_id, doiless_publication_id, None, None))
                continue

            url = f"https://doi.org/{doiless_publication_id}"
            try:
                response = requests.get(url, timeout=REQUEST_TIMEOUT)
                checks.append(DoiCheck(citation, publication_id, doiless_publication_id, response.status_code, None))
            except requests.RequestException as exc:
                checks.append(DoiCheck(citation, publication_id, doiless_publication_id, None, str(exc)))
        return checks

    return cached_lint_result(tool_source, "doi", calculate)


class DoiEmpty(Linter):
    @classmethod
    def lint(cls, tool_source: "ToolSource", lint_ctx: "LintContext"):
        for check in _doi_checks(tool_source):
            if not check.doiless_publication_id:
                lint_ctx.error("Empty DOI citation", linter=cls.name(), node=check.citation)


class DoiInvalid(Linter):
    @classmethod
    def lint(cls, tool_source: "ToolSource", lint_ctx: "LintContext"):
        for check in _doi_checks(tool_source):
            if check.status_code == 404:
                lint_ctx.error(
                    f"{check.publication_id} is not a valid DOI",
                    linter=cls.name(),
                    node=check.citation,
                )


class DoiPrefix(Linter):
    @classmethod
    def lint(cls, tool_source: "ToolSource", lint_ctx: "LintContext"):
        for check in _doi_checks(tool_source):
            if check.status_code == 200 and check.publication_id != check.doiless_publication_id:
                lint_ctx.error(
                    f"{check.publication_id} is valid, but Galaxy expects DOI without 'doi:' prefix",
                    linter=cls.name(),
                    node=check.citation,
                )


class DoiUnexpectedResponse(Linter):
    @classmethod
    def lint(cls, tool_source: "ToolSource", lint_ctx: "LintContext"):
        for check in _doi_checks(tool_source):
            if check.error:
                url = f"https://doi.org/{check.doiless_publication_id}"
                lint_ctx.warn(f"Error '{check.error}' accessing {url}", linter=cls.name(), node=check.citation)
            elif check.status_code is not None and check.status_code not in [200, 404]:
                lint_ctx.warn(
                    f"doi.org returned unexpected status code {check.status_code}",
                    linter=cls.name(),
                    node=check.citation,
                )


class DoiValid(Linter):
    @classmethod
    def lint(cls, tool_source: "ToolSource", lint_ctx: "LintContext"):
        for check in _doi_checks(tool_source):
            if check.status_code == 200 and check.publication_id == check.doiless_publication_id:
                lint_ctx.info(
                    f"{check.publication_id} is a valid DOI",
                    linter=cls.name(),
                    node=check.citation,
                )
