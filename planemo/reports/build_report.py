import base64

from galaxy.util import strip_control_characters
from jinja2 import (
    Environment,
    PackageLoader,
)
from pkg_resources import resource_string

TITLE = "Results (powered by Planemo)"


def build_report(structured_data, report_type="html", execution_type="Test", **kwds):
    """Use report_{report_type}.tpl to build page for report."""
    environment = dict(
        title=TITLE,
        raw_data=structured_data,
    )

    __fix_test_ids(environment)
    environment = __inject_summary(environment)

    if report_type == "html":
        # The HTML report format needs a lot of extra, custom data.
        # IMO, this seems to suggest it should be embedded.
        environment["title"] = None
        environment["execution_type"] = execution_type
        markdown = template_data(environment, "report_markdown.tpl")
        environment["title"] = " ".join((environment["execution_type"], TITLE))
        environment["raw_data"] = base64.b64encode(markdown.encode("utf-8")).decode("utf-8")
        environment.update(
            {
                "custom_style": __style("custom.css"),
                "custom_script": __script("custom"),
                "bootstrap_style": __style("bootstrap.min.css"),
                "jquery_script": __script("jquery.min"),
                "bootstrap_script": __script("bootstrap.min"),
                "markdown_it_script": __script("markdown-it.min"),
            }
        )

    return template_data(environment, "report_%s.tpl" % report_type)


def template_data(environment, template_name, **kwds):
    """Build an arbitrary templated page."""
    env_kwargs = {}
    if template_name == "report_markdown.tpl":
        env_kwargs["keep_trailing_newline"] = True
        env_kwargs["trim_blocks"] = True
    env = Environment(loader=PackageLoader("planemo", "reports"), **env_kwargs)
    env.filters["strip_control_characters"] = lambda x: strip_control_characters(x) if x else x
    template = env.get_template(template_name)
    return template.render(**environment)


def __fix_test_ids(environment):
    for test in environment["raw_data"]["tests"]:
        test_data = test.get("data")
        if test_data and test_data.get("tool_id"):
            test["id"] = f"{test_data['tool_id']} (Test #{test_data['test_index'] + 1})"


def __inject_summary(environment):
    total = 0
    errors = 0
    failures = 0
    skips = 0
    for execution in environment["raw_data"]["tests"]:
        total += 1
        test_data = execution.get("data")
        if test_data:
            status = test_data.get("status")
            if status == "error":
                errors += 1
            elif status == "failure":
                failures += 1
            elif status == "skipped":
                skips += 1
    environment["raw_data"]["results"] = {
        "total": total,
        "errors": errors,
        "failures": failures,
        "skips": skips,
    }
    if "suitename" not in environment:
        environment["raw_data"]["suitename"] = TITLE
    return environment


def __style(filename):
    resource = __load_resource(filename)
    return "<style>%s</style>" % resource


def __script(short_name):
    resource = __load_resource("%s.js" % short_name)
    return "<script>%s</script>" % resource


def __load_resource(name):
    return resource_string(__name__, name).decode("UTF-8")
