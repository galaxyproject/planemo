import os

from galaxy.util import unicodify

from planemo import shed_lint
from planemo.xml import validation
from .test_utils import (
    skip_unless_executable,
    skip_unless_module,
    TEST_DIR,
)


@skip_unless_module("lxml")
def test_lxml_validation():
    lxml_xsd_validator = validation.LxmlValidator()
    _check_validator(lxml_xsd_validator)


@skip_unless_executable("xmllint")
def test_xmllint_validation():
    xmllint_xsd_validator = validation.XmllintValidator()
    _check_validator(xmllint_xsd_validator)


def test_tool_dependencies_validation():
    _assert_validates(shed_lint.TOOL_DEPENDENCIES_XSD, _path("tool_dependencies_good_1.xml"))
    _assert_validates(shed_lint.TOOL_DEPENDENCIES_XSD, _path("tool_dependencies_good_2.xml"))


def test_repository_dependencies_validation():
    _assert_validates(shed_lint.REPO_DEPENDENCIES_XSD, _path("repository_dependencies.xml"))


def _check_validator(xsd_validator):
    _assert_validates(
        _path("xsd_schema_1.xsd"),
        _path("xml_good_1.xml"),
        xsd_validator,
    )
    _assert_validates(
        _path("xsd_schema_1.xsd"),
        _path("xml_good_2.xml"),
        xsd_validator,
    )

    result = xsd_validator.validate(_path("xsd_schema_1.xsd"), _path("xml_bad_1.xml"))
    assert not result.passed
    output = unicodify(result.output)
    assert "not_command" in output, output


def _assert_validates(schema, target, xsd_validator=None):
    if xsd_validator is None:
        xsd_validator = validation.get_validator()
    result = xsd_validator.validate(schema, target)
    assert result.passed, result.output


def _path(filename):
    # TODO: move all test data into tests/data
    return os.path.join(TEST_DIR, filename)
