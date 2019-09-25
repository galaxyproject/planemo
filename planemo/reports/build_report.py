import json

from galaxy.util import strip_control_characters
from jinja2 import Environment, PackageLoader
from pkg_resources import resource_string

TITLE = "Tool Test Results (powered by Planemo)"
env = Environment(loader=PackageLoader('planemo', 'reports'))
env.filters['strip_control_characters'] = lambda x: strip_control_characters(x) if x else x


def build_report(structured_data, report_type="html", **kwds):
    """ Use report_{report_type}.tpl to build page for report.
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
            'json': json,
        })

    environment = __inject_summary(environment)

    return template_data(environment, 'report_%s.tpl' % report_type)


def template_data(environment, template_name="report_html.tpl", **kwds):
    """Build an arbitrary templated page.
    """
    template = env.get_template(template_name)
    return template.render(**environment)


def __inject_summary(environment):
    if 'results' not in environment['raw_data']:
        total = 0
        errors = 0
        failures = 0
        skips = 0
        for test in environment['raw_data']['tests']:
            total += 1
            test_data = test.get('data')
            if test_data:
                status = test_data.get('status')
                if status != 'ok':
                    errors += 1
        environment['raw_data']['results'] = {
            'total': total,
            'errors': errors,
            'failures': failures,
            'skips': skips,
        }
    if 'suitename' not in environment:
        environment['raw_data']['suitename'] = TITLE
    return environment


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
