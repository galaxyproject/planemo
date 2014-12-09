import json
from string import Template
from pkg_resources import resource_string

TITLE = "Tool Test Results (powered by Planemo)"
LINKS_HTML = """
    <li><a href="https://galaxyproject.org">Galaxy</a></li>
    <li><a href="https://planemo.readthedocs.com">Planemo</a></li>
"""


def build_report(structured_data, **kwds):
    """ Use report_template.html to build HTML page for report.
    """
    custom_style = __style("custom.css")
    custom_script = __script("custom")
    bootstrap_style = __style("bootstrap.min.css")
    jquery_script = __script("jquery.min")
    bootstrap_script = __script("bootstrap.min")

    environment = dict(
        custom_style=custom_style,
        custom_script=custom_script,
        bootstrap_style=bootstrap_style,
        jquery_script=jquery_script,
        bootstrap_script=bootstrap_script,
        title=TITLE,
        links=LINKS_HTML,
        json_test_data=json.dumps(structured_data),
    )
    template = Template(__load_resource("report_template.html"))
    return template.safe_substitute(environment)


def __style(filename):
    resource = __load_resource(filename)
    return "<style>%s</style>" % resource


def __script(short_name):
    resource = __load_resource("%s.js" % short_name)
    return "<script>%s</script>" % resource


def __load_resource(name):
    return resource_string(
        __name__, name
    ).decode('UTF-8')
