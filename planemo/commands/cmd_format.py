"""Module describing the planemo ``format`` command."""

import difflib
import os

import click
from lxml import etree

from planemo import options
from planemo.cli import (
    command_function,
    PlanemoCliContext,
)
from planemo.io import (
    error,
    info,
    warn,
)


def format_xml(content: str, tab_size: int = 4) -> str:
    """Format XML content with consistent indentation."""
    parser = etree.XMLParser(strip_cdata=False)
    xml = etree.fromstring(content, parser=parser)
    etree.indent(xml, space=" " * tab_size)
    return etree.tostring(xml, pretty_print=True, encoding=str)


def find_xml_files(paths, recursive=False):
    """Find all .xml files in the given paths."""
    xml_files = []
    for path in paths:
        if os.path.isfile(path):
            if path.endswith(".xml"):
                xml_files.append(path)
        elif os.path.isdir(path):
            if recursive:
                for root, _dirs, files in os.walk(path):
                    for f in sorted(files):
                        if f.endswith(".xml"):
                            xml_files.append(os.path.join(root, f))
            else:
                for f in sorted(os.listdir(path)):
                    if f.endswith(".xml"):
                        xml_files.append(os.path.join(path, f))
    return xml_files


@click.command("format")
@options.optional_tools_arg(multiple=True)
@options.recursive_option()
@click.option(
    "--dry-run",
    is_flag=True,
    default=False,
    help="Show diff of changes without modifying files.",
)
@click.option(
    "--tab-size",
    type=int,
    default=4,
    help="Number of spaces per indentation level (default: 4).",
)
@command_function
def cli(ctx: PlanemoCliContext, paths, **kwds):
    """Format XML files with consistent indentation."""
    if not paths:
        paths = [os.getcwd()]

    xml_files = find_xml_files(paths, recursive=kwds["recursive"])
    if not xml_files:
        warn("No XML files found.")
        ctx.exit(0)

    dry_run = kwds["dry_run"]
    tab_size = kwds["tab_size"]
    changed = 0

    for xml_path in xml_files:
        try:
            with open(xml_path) as f:
                original = f.read()
        except Exception as e:
            error("Could not read %s: %s", xml_path, e)
            continue

        try:
            formatted = format_xml(original, tab_size=tab_size)
        except etree.XMLSyntaxError as e:
            warn("Skipping %s (XML syntax error: %s)", xml_path, e)
            continue

        if formatted != original:
            changed += 1
            if dry_run:
                diff = difflib.unified_diff(
                    original.splitlines(keepends=True),
                    formatted.splitlines(keepends=True),
                    fromfile=xml_path,
                    tofile=xml_path,
                )
                click.echo("".join(diff))
            else:
                with open(xml_path, "w") as f:
                    f.write(formatted)
                info("Formatted %s", xml_path)
        else:
            info("%s already formatted", xml_path)

    if dry_run:
        info("%d file(s) would be reformatted.", changed)
    else:
        info("%d file(s) reformatted.", changed)
    ctx.exit(0)
