TAG_ORDER = [
    'description',
    'macros',
    'requirements',
    'code',
    'stdio',
    'version_command',
    'command',
    'inputs',
    'configfiles',
    'outputs',
    'tests',
    'help',
    'citations',
]


# Ensure the XML blocks appear in the correct order prescribed
# by the tool author best practices.
def lint_xml_ordering(tool_xml, lint_ctx):
    last_tag = None
    last_key = None
    for elem in list(tool_xml.getroot()):
        tag = elem.tag
        if tag in TAG_ORDER:
            key = TAG_ORDER.index(tag)
            if last_key:
                if last_key > key:
                    lint_ctx.warn("Best practice violation [%s] elements should come before [%s]" % (tag, last_tag))
            last_tag = tag
            last_key = key
        else:
            lint_ctx.info("Unknown tag [%s] encoutered, this may result in a warning in the future." % tag)
