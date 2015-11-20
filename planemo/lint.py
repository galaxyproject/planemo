import os

from six.moves.urllib.request import (
    urlopen,
)
from six.moves.urllib.error import (
    HTTPError,
    URLError,
)
from planemo.xml import validation
from planemo.shed import find_urls_for_xml


def lint_xsd(lint_ctx, schema_path, path):
    name = os.path.basename(path)
    validator = validation.get_validator(require=True)
    validation_result = validator.validate(schema_path, path)
    if not validation_result.passed:
        msg = "Invalid %s found. Errors [%s]"
        msg = msg % (name, validation_result.output)
        lint_ctx.error(msg)
    else:
        lint_ctx.info("%s found and appears to be valid XML" % name)


def lint_urls(root, lint_ctx):
    urls = find_urls_for_xml(root)

    def validate_url(url, lint_ctx):
        try:
            handle = urlopen(url)
            handle.read(100)
            lint_ctx.info("URL OK %s" % url)
        except HTTPError as e:
            lint_ctx.error("HTTP Error %s accessing %s" % (e.code, url))
        except URLError as e:
            lint_ctx.error("URL Error %s accessing %s" % (str(e), url))

    for url in urls:
        validate_url(url, lint_ctx)
