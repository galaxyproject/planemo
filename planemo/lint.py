"""Utilities to help linting various targets."""
from __future__ import absolute_import

import os

import requests

from galaxy.tools.lint import LintContext
from six.moves.urllib.error import (
    HTTPError,
    URLError,
)
from six.moves.urllib.request import (
    Request,
    urlopen,
)

import planemo.linters.biocontainer_registered
import planemo.linters.conda_requirements
import planemo.linters.doi
import planemo.linters.urls
import planemo.linters.xsd

from planemo.io import error
from planemo.shed import find_urls_for_xml
from planemo.xml import validation


def build_lint_args(ctx, **kwds):
    """Handle common report, error, and skip linting arguments."""
    report_level = kwds.get("report_level", "all")
    fail_level = kwds.get("fail_level", "warn")
    skip = kwds.get("skip", None)
    if skip is None:
        skip = ctx.global_config.get("lint_skip", "")
        if isinstance(skip, list):
            skip = ",".join(skip)

    skip_types = [s.strip() for s in skip.split(",")]
    lint_args = dict(
        level=report_level,
        fail_level=fail_level,
        extra_modules=_lint_extra_modules(**kwds),
        skip_types=skip_types,
    )
    return lint_args


# TODO: Move this back to tool_lint.
def _lint_extra_modules(**kwds):
    linters = []
    if kwds.get("xsd", True):
        linters.append(planemo.linters.xsd)

    if kwds.get("doi", False):
        linters.append(planemo.linters.doi)

    if kwds.get("urls", False):
        linters.append(planemo.linters.urls)

    if kwds.get("conda_requirements", False):
        linters.append(planemo.linters.conda_requirements)

    if kwds.get("biocontainer", False):
        linters.append(planemo.linters.biocontainer_registered)

    return linters


def setup_lint(ctx, **kwds):
    """Setup lint_args and lint_ctx to begin linting a target."""
    lint_args = build_lint_args(ctx, **kwds)
    lint_ctx = LintContext(lint_args["level"])
    return lint_args, lint_ctx


def handle_lint_complete(lint_ctx, lint_args, failed=False):
    """Complete linting of a target and decide exit code."""
    if not failed:
        failed = lint_ctx.failed(lint_args["fail_level"])
    if failed:
        error("Failed linting")
    return 1 if failed else 0


def lint_dois(tool_xml, lint_ctx):
    """Find referenced DOIs and check they have valid with http://dx.doi.org."""
    dois = find_dois_for_xml(tool_xml)
    for publication in dois:
        is_doi(publication, lint_ctx)


def find_dois_for_xml(tool_xml):
    dois = []
    for element in tool_xml.getroot().findall("citations"):
        for citation in list(element):
            if citation.tag == 'citation' and citation.attrib.get('type', '') == 'doi':
                dois.append(citation.text)
    return dois


def is_doi(publication_id, lint_ctx):
    """Check if dx.doi knows about the ``publication_id``."""
    base_url = "http://dx.doi.org"
    doiless_publication_id = publication_id.split("doi:", 1)[-1]
    url = "%s/%s" % (base_url, doiless_publication_id)
    r = requests.get(url)
    if r.status_code == 200:
        if publication_id != doiless_publication_id:
            lint_ctx.error("%s is valid, but Galaxy expects DOI without 'doi:' prefix" % publication_id)
        else:
            lint_ctx.info("%s is a valid DOI" % publication_id)
    elif r.status_code == 404:
        lint_ctx.error("%s is not a valid DOI" % publication_id)
    else:
        lint_ctx.warn("dx.doi returned unexpected status code %d" % r.status_code)


def lint_xsd(lint_ctx, schema_path, path):
    """Lint XML at specified path with supplied schema."""
    name = os.path.basename(path)
    validator = validation.get_validator(require=True)
    validation_result = validator.validate(schema_path, path)
    if not validation_result.passed:
        msg = "Invalid %s found. Errors [%s]"
        msg = msg % (name, validation_result.output)
        lint_ctx.error(msg)
    else:
        lint_ctx.info("File validates against XML schema.")


def lint_urls(root, lint_ctx):
    """Find referenced URLs and verify they are valid."""
    urls, docs = find_urls_for_xml(root)

    # This is from Google Chome 53.0.2785.143, current at time of writing:
    BROWSER_USER_AGENT = "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_11_6) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/53.0.2785.143 Safari/537.36"

    def validate_url(url, lint_ctx, user_agent=None):
        is_valid = True
        if user_agent:
            req = Request(url, headers={"User-Agent": user_agent})
        else:
            req = url
        try:
            handle = urlopen(req)
            handle.read(100)
        except HTTPError as e:
            if e.code == 429:
                # too many requests
                pass
            else:
                is_valid = False
                lint_ctx.error("HTTP Error %s accessing %s" % (e.code, url))
        except URLError as e:
            is_valid = False
            lint_ctx.error("URL Error %s accessing %s" % (str(e), url))
        if is_valid:
            lint_ctx.info("URL OK %s" % url)

    for url in urls:
        validate_url(url, lint_ctx)
    for url in docs:
        validate_url(url, lint_ctx, BROWSER_USER_AGENT)


__all__ = (
    "build_lint_args",
    "handle_lint_complete",
    "lint_dois",
    "lint_urls",
    "lint_xsd",
)
