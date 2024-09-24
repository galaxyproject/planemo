""" Tool linting module that lints Galaxy tools for their URLs
"""

from typing import TYPE_CHECKING
from urllib.request import urlopen

import requests
from galaxy.tool_util.lint import Linter

from planemo.shed import _find_urls_in_text

if TYPE_CHECKING:
    from galaxy.tool_util.lint import LintContext
    from galaxy.tool_util.parser.interface import ToolSource
    from galaxy.util import ElementTree

BROWSER_USER_AGENT = "Mozilla/5.0 (Macintosh; Intel Mac OS X 11_1) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/87.0.4280.141 Safari/537.36"


def find_urls_in_help(root: "ElementTree"):
    for help in root.findall("help"):
        for url in _find_urls_in_text(help.text):
            yield url[0], help


class URLInaccessibleHttp(Linter):
    @classmethod
    def lint(cls, tool_source: "ToolSource", lint_ctx: "LintContext"):
        tool_xml = getattr(tool_source, "xml_tree", None)
        if not tool_xml:
            return
        for url, help in find_urls_in_help(tool_xml):
            if url.startswith("http://") or url.startswith("https://"):
                headers = {"User-Agent": BROWSER_USER_AGENT, "Accept": "*/*"}
                r = None
                try:
                    r = requests.get(url, headers=headers, stream=True)
                    r.raise_for_status()
                    next(r.iter_content(1000))
                except Exception as e:
                    if r is not None and r.status_code == 429:
                        # too many requests
                        pass
                    if r is not None and r.status_code in [403, 503] and "cloudflare" in r.text:
                        # CloudFlare protection block
                        pass
                    else:
                        lint_ctx.error(f"Error '{e}' accessing {url}", linter=cls.name(), node=help)


class URLInaccessible(Linter):
    @classmethod
    def lint(cls, tool_source: "ToolSource", lint_ctx: "LintContext"):
        tool_xml = getattr(tool_source, "xml_tree", None)
        if not tool_xml:
            return
        for url, help in find_urls_in_help(tool_xml):
            if not url.startswith("http://") and not url.startswith("https://"):
                try:
                    with urlopen(url) as handle:
                        handle.read(100)
                except Exception as e:
                    lint_ctx.error(f"Error '{e}' accessing {url}", linter=cls.name(), node=help)


class URLValid(Linter):
    @classmethod
    def lint(cls, tool_source: "ToolSource", lint_ctx: "LintContext"):
        tool_xml = getattr(tool_source, "xml_tree", None)
        if not tool_xml:
            return
        for url, help in find_urls_in_help(tool_xml):
            is_valid = True
            if url.startswith("http://") or url.startswith("https://"):
                headers = {"User-Agent": BROWSER_USER_AGENT, "Accept": "*/*"}
                r = None
                try:
                    r = requests.get(url, headers=headers, stream=True)
                    r.raise_for_status()
                    next(r.iter_content(1000))
                except Exception:
                    if r is None or (
                        r.status_code != 429 and not (r.status_code in [403, 503] and "cloudflare" in r.text)
                    ):
                        is_valid = False
            else:
                try:
                    with urlopen(url) as handle:
                        handle.read(100)
                except Exception:
                    is_valid = False
            if is_valid:
                lint_ctx.info("URL OK %s" % url, linter=cls.name(), node=help)
