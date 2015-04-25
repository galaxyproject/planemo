try:
    from jinja2 import Template
except ImportError:
    Template = None

NO_JINJA2_MESSAGE = ("This functionality requires Jinja2 but this library is "
                     "unavailable. Install with `pip install jinja2`.")


def render(template_str, **kwds):
    """ Use jinja2 to render a template
    """
    if Template is None:
        raise Exception(NO_JINJA2_MESSAGE)
    template = Template(template_str)
    contents = template.render(**kwds)
    return contents
