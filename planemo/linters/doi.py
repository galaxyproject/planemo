""" Tool linting module that lints Galaxy tools for their DOIs (if a DOI type citation is present)
"""

from typing import TYPE_CHECKING

import requests
from galaxy.tool_util.lint import Linter

if TYPE_CHECKING:
    from galaxy.tool_util.lint import LintContext
    from galaxy.tool_util.parser.interface import ToolSource
    from galaxy.util import ElementTree


class DoiEmptyNone(Linter):
    @classmethod
    def lint(cls, tool_source: "ToolSource", lint_ctx: "LintContext"):
        tool_xml = getattr(tool_source, "xml_tree", None)
        if not tool_xml:
            return
        for citation, *_ in _doi_citations(tool_xml):
            if citation.text is None:
                lint_ctx.error("Empty DOI citation", linter=cls.name(), node=citation)


class DoiEmpty(Linter):
    @classmethod
    def lint(cls, tool_source: "ToolSource", lint_ctx: "LintContext"):
        tool_xml = getattr(tool_source, "xml_tree", None)
        if not tool_xml:
            return
        for citation, publication_id, doiless_publication_id in _doi_citations(tool_xml):
            if citation.text is None:
                continue
            if not doiless_publication_id:
                lint_ctx.error("Empty DOI citation", linter=cls.name(), node=citation)


class DoiValid(Linter):
    @classmethod
    def lint(cls, tool_source: "ToolSource", lint_ctx: "LintContext"):
        tool_xml = getattr(tool_source, "xml_tree", None)
        if not tool_xml:
            return
        for citation, publication_id, doiless_publication_id in _doi_citations(tool_xml):
            if citation.text is None or not doiless_publication_id:
                continue
            url = f"https://doi.org/{doiless_publication_id}"
            r = requests.get(url)
            if r.status_code == 200 and publication_id == doiless_publication_id:
                lint_ctx.info("%s is a valid DOI" % publication_id, linter=cls.name(), node=citation)


class DoiValidWithDoi(Linter):
    @classmethod
    def lint(cls, tool_source: "ToolSource", lint_ctx: "LintContext"):
        tool_xml = getattr(tool_source, "xml_tree", None)
        if not tool_xml:
            return
        for citation, publication_id, doiless_publication_id in _doi_citations(tool_xml):
            if citation.text is None or not doiless_publication_id:
                continue
            url = f"https://doi.org/{doiless_publication_id}"
            r = requests.get(url)
            if r.status_code == 200 and publication_id != doiless_publication_id:
                lint_ctx.error(
                    "%s is valid, but Galaxy expects DOI without 'doi:' prefix" % publication_id,
                    linter=cls.name(),
                    node=citation,
                )


class DoiInvalid(Linter):
    @classmethod
    def lint(cls, tool_source: "ToolSource", lint_ctx: "LintContext"):
        tool_xml = getattr(tool_source, "xml_tree", None)
        if not tool_xml:
            return
        for citation, publication_id, doiless_publication_id in _doi_citations(tool_xml):
            if citation.text is None or not doiless_publication_id:
                continue
            url = f"https://doi.org/{doiless_publication_id}"
            r = requests.get(url)
            if r.status_code == 404:
                lint_ctx.error("%s is not a valid DOI" % publication_id, linter=cls.name(), node=citation)


class DoiUnexpectedResponse(Linter):
    @classmethod
    def lint(cls, tool_source: "ToolSource", lint_ctx: "LintContext"):
        tool_xml = getattr(tool_source, "xml_tree", None)
        if not tool_xml:
            return
        for citation, publication_id, doiless_publication_id in _doi_citations(tool_xml):
            if citation.text is None or not doiless_publication_id:
                continue
            url = f"https://doi.org/{doiless_publication_id}"
            r = requests.get(url)
            if r.status_code not in [200, 400]:
                lint_ctx.warn(
                    "dx.doi returned unexpected status code %d" % r.status_code, linter=cls.name(), node=citation
                )


def _doi_citations(tool_xml: "ElementTree"):
    for element in tool_xml.getroot().findall("citations"):
        for citation in list(element):
            if citation.tag == "citation" and citation.attrib.get("type", "") == "doi":
                if citation.text is None:
                    publication_id = doiless_publication_id = None
                else:
                    publication_id = citation.text.strip()
                    doiless_publication_id = publication_id.split("doi:", 1)[-1]
                yield citation, publication_id, doiless_publication_id
