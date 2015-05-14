from __future__ import print_function

import os
import sys
from xml.etree import ElementTree

from planemo.xml import diff


def diff_and_remove(working, label_a, label_b, f):
    assert label_a != label_b
    special = ["tool_dependencies.xml", "repository_dependencies.xml"]
    deps_diff = 0
    # Could walk either A or B; will only compare if in same relative location
    for dirpath, dirnames, filenames in os.walk(os.path.join(working, label_a)):
        for filename in filenames:
            if filename in special:
                a = os.path.join(dirpath, filename)
                b = os.path.join(working, label_b, os.path.relpath(a, os.path.join(working, label_a)))
                if os.path.exists(a) and os.path.exists(b):
                    deps_diff &= _shed_diff(a, b, f)
                    os.remove(a)
                    os.remove(b)
    return deps_diff


def _shed_diff(file_a, file_b, f=sys.stdout):
    xml_a = ElementTree.parse(file_a).getroot()
    xml_b = ElementTree.parse(file_b).getroot()
    _strip_shed_attributes(xml_a)
    _strip_shed_attributes(xml_b)
    return diff.diff(xml_a, xml_b, reporter=f.write)


def _strip_shed_attributes(xml_element):
    if xml_element.tag == "repository":
        _remove_attribs(xml_element)
    children = xml_element.getchildren()
    if len(children) > 0:
        for child in children:
            _strip_shed_attributes(child)


def _remove_attribs(xml_element):
    for attrib in ["changeset_revision", "toolshed"]:
        if attrib in xml_element.attrib:
            del xml_element.attrib[attrib]
