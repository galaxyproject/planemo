from __future__ import print_function

import os
import sys
from xml.etree import ElementTree

from planemo.xml import diff


def diff_and_remove(working, label_a, label_b, f):
    a_deps = os.path.join(working, label_a, "tool_dependencies.xml")
    b_deps = os.path.join(working, label_b, "tool_dependencies.xml")
    a_repos = os.path.join(working, label_a, "repository_dependencies.xml")
    b_repos = os.path.join(working, label_b, "repository_dependencies.xml")

    deps_diff = 0
    if os.path.exists(a_deps) and os.path.exists(b_deps):
        deps_diff = _shed_diff(a_deps, b_deps, f)
        os.remove(a_deps)
        os.remove(b_deps)

    repos_diff = 0
    if os.path.exists(a_repos) and os.path.exists(b_repos):
        repos_diff = _shed_diff(a_repos, b_repos, f)
        os.remove(a_repos)
        os.remove(b_repos)

    return deps_diff and repos_diff


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
