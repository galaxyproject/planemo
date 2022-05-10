import os
import sys
from xml.etree import ElementTree

from planemo.xml.diff import diff
from .test_utils import TEST_DIR


def test_diff():
    # Check with and without reporter.
    _check(sys.stdout.write)
    _check(None)


def _check(reporter):
    assert not diff(ElementTree.fromstring("<moo>cow</moo>"), ElementTree.fromstring("<moo>cow</moo>"), reporter)
    assert diff(ElementTree.fromstring("<moo>cow</moo>"), ElementTree.fromstring("<moo>cow1</moo>"), reporter)
    assert diff(ElementTree.fromstring("<moo><c/><c1/></moo>"), ElementTree.fromstring("<moo><c/></moo>"), reporter)
    assert diff(ElementTree.fromstring("<moo><c/></moo>"), ElementTree.fromstring("<moo><c/>a</moo>"), reporter)
    assert diff(ElementTree.fromstring("<moo>*</moo>"), ElementTree.fromstring("<moo>*a</moo>"), reporter)
    assert not diff(_root("repository_dependencies.xml"), _root("repository_dependencies.xml"), reporter)
    assert not diff(_root("repository_dependencies_shed.xml"), _root("repository_dependencies_shed.xml"), reporter)
    assert diff(_root("repository_dependencies.xml"), _root("repository_dependencies_shed.xml"), reporter)
    assert diff(_root("repository_dependencies.xml"), _root("tool_dependencies_good_1.xml"), reporter)
    assert diff(_root("tool_dependencies_good_1.xml"), _root("tool_dependencies_good_2.xml"), reporter)


def _root(*args):
    return ElementTree.parse(os.path.join(TEST_DIR, *args)).getroot()
