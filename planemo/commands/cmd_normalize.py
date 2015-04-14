from xml.etree import ElementTree

import click

from planemo.cli import pass_context
from planemo import options

from galaxy.tools.loader import (
    load_tool,
    raw_tool_xml_tree,
)
from galaxy.tools.linters.xml_order import TAG_ORDER


@click.command('normalize')
@options.required_tool_arg()
@click.option(
    "--expand_macros",
    is_flag=True,
    help=("Expand macros while normalizing tool XML - useful to see how "
          "macros are evaluated.")
)
@click.option(
    "--skip_reorder",
    is_flag=True,
    help=("Planemo will reorder top-level tool blocks according to tool "
          "development best practices as part of this command, this flag "
          "will disable that behavior.")
)
@click.option(
    "--skip_reindent",
    is_flag=True,
    help=("Planemo will reindent the XML according to tool development "
          "best practices as part of this command, this flag will disable "
          "that behavior.")
)
@pass_context
def cli(ctx, path, expand_macros=False, **kwds):
    """Generate normalized tool XML from input (breaks formatting).

    This will break the formatting of your tool and is currently only intended
    for viewing macro expansions for for use with XSD validation (see
    https://github.com/JeanFred/Galaxy-XSD for instance). Please do not use
    the output as is - it frequently makes tool less readable not more.

    The top-level blocks will be reordered and whitespace fixed according to
    the tool development best practices outlined on the Galaxy wiki.

    ::

        % # Print normalized version of tool.
        % planemo normalize tool.xml
        <tool>
        ...
        % # Print a variant of tool with all macros expanded out, useful for
        % # debugging complex macros.
        % planemo normalize --expand_macros tool.xml
        <tool>
        ...
    """
    if expand_macros:
        tree = load_tool(path)
    else:
        tree = raw_tool_xml_tree(path)

    root = tree.getroot()
    if not kwds.get("skip_reorder", False):
        last_index = len(TAG_ORDER)
        data = []
        for elem in root:
            tag = elem.tag
            if tag in TAG_ORDER:
                key = TAG_ORDER.index(tag)
            else:
                key = last_index
                last_index += 1
            data.append((key, elem))
        data.sort()
        root[:] = [item[-1] for item in data]
    if not kwds.get("skip_reindent", False):
        _indent(root)
    ElementTree.dump(root)


def _indent(elem, level=0):
    # http://stackoverflow.com/questions/749796/pretty-printing-xml-in-python
    i = "\n" + level * "    "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "    "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            _indent(elem, level + 1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i
