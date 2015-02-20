""" Tool linting module that lints Galaxy tool against experimental XSD.
"""
import copy
import os
import tempfile

from planemo.xml import XSDS_PATH
import planemo.lint

TOOL_XSD = os.path.join(XSDS_PATH, "tool", "galaxy.xsd")


def lint_tool_xsd(root, lint_ctx):
    """ Write a temp file out and lint it.
    """
    with tempfile.NamedTemporaryFile() as tf:
        _clean_root(root).write(tf.name)
        planemo.lint.lint_xsd(lint_ctx, TOOL_XSD, tf.name)


def _clean_root(root):
    """ XSD assumes macros have been expanded, so remove them.
    """
    clean_root = copy.deepcopy(root)
    to_remove = []
    for macros_el in clean_root.findall("macros"):
        to_remove.append(macros_el)
    for macros_el in to_remove:
        clean_root.getroot().remove(macros_el)
    return clean_root
