import json
from pkg_resources import resource_string
from jinja2 import Environment, PackageLoader
env = Environment(loader=PackageLoader('planemo', 'reports'))
TITLE = "Tool Test Results (powered by Planemo)"


def build_report(structured_data, report_type="html", **kwds):
    """ Use report_template.html to build HTML page for report.
    """
    environment = dict(
        title=TITLE,
        raw_data=structured_data,
    )

    if report_type == 'html':
        # The HTML report format needs a lot of extra, custom data.
        # IMO, this seems to suggest it should be embedded.
        environment.update({
            'custom_style': __style("custom.css"),
            'custom_script': __script("custom"),
            'bootstrap_style': __style("bootstrap.min.css"),
            'jquery_script': __script("jquery.min"),
            'bootstrap_script': __script("bootstrap.min"),
            'json_test_data': json.dumps(structured_data),
        })

    template = env.get_template('report_%s.tpl' % report_type)
    return template.render(**environment)


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
