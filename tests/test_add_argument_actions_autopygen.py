import ast
import os
import sys
from typing import (
    Dict,
    Optional,
    Set,
    Tuple,
)

import pytest

from planemo.autopygen.argument_parser_conversion import (
    command_from_decoy,
    obtain_parser,
    xml_from_decoy,
    xml_to_string,
)
from planemo.autopygen.param_info import ParamInfo
from tests.test_utils import (
    assert_equal,
    load_function_body,
    TEST_AUTOPYGEN_DATA,
)


def test_no_action_positional():
    inputs, _, _ = extract_xml("none_positional")
    assert_equal(inputs, '<param argument="test" type="text" optional="false" label="test"/>\n')


def test_no_action():
    inputs, _, _ = extract_xml("none")
    assert_equal(inputs, '<param argument="--test" type="text" optional="true" label="test"/>\n')


def test_store():
    inputs, _, _ = extract_xml("store")
    assert_equal(inputs, '<param argument="--test" type="text" optional="true" label="test"/>\n')


def test_store_with_default():
    inputs, _, _ = extract_xml("store_with_default")
    assert_equal(inputs, '<param argument="--test" type="text" value="foo" optional="true" label="test"/>\n')


STORE_WITH_DEFAULT_CHOICES = (
    '<param argument="--test" type="select" optional="true" label="test">\n'
    '    <option value="foo">Foo</option>\n'
    '    <option value="bar" selected="true">Bar</option>\n'
    '    <option value="goo">Goo</option>\n'
    "</param>\n"
)


def test_store_with_default_choices():
    inputs, _, _ = extract_xml("store_with_default_choices")
    assert_equal(inputs, STORE_WITH_DEFAULT_CHOICES)


def test_store_const():
    inputs, _, _ = extract_xml("store_const")
    assert_equal(
        inputs,
        '<param argument="--test" type="boolean" truevalue="--test" falsevalue="" checked="false" ' 'label="test"/>\n',
    )


def test_store_true():
    inputs, _, _ = extract_xml("store_true")
    assert_equal(
        inputs,
        '<param argument="--test" type="boolean" truevalue="--test" falsevalue="" checked="false" ' 'label="test"/>\n',
    )


def test_append():
    inputs, _, _ = extract_xml("append")
    expected = (
        '<repeat name="test_repeat" title="test_repeat">\n'
        '    <param argument="--test" type="text" optional="true" label="test"/>\n'
        "</repeat>\n"
    )

    assert_equal(inputs, expected)


def test_append_const():
    inputs, _, _ = extract_xml("append_const")
    expected = (
        '<repeat name="test_repeat" title="test_repeat">\n'
        '    <param argument="--test" type="boolean" truevalue="--test" falsevalue="" checked="false" '
        'label="test"/>\n'
        "</repeat>\n"
    )
    assert_equal(inputs, expected)


def test_count():
    inputs, _, _ = extract_xml("count")
    expected = (
        '<repeat name="test_repeat" title="test_repeat">\n'
        '    <param argument="--test" type="boolean" truevalue="--test" falsevalue="" checked="false" '
        'label="test"/>\n'
        "</repeat>\n"
    )
    assert_equal(inputs, expected)


def test_version():
    inputs, _, param_info = extract_xml("version")

    assert_equal(inputs, "")
    assert_equal(type(param_info), ParamInfo)
    assert_equal(param_info.param_type.is_version, True)


@pytest.mark.skipif(sys.version_info < (3, 8), reason="Extend command is supported from python 3.8")
def test_extend():
    inputs, _, _ = extract_xml("extend")
    expected = (
        '<repeat name="test_repeat" title="test_repeat">\n'
        '    <param argument="--test" type="text" optional="true" label="test"/>\n'
        "</repeat>\n"
    )
    assert_equal(inputs, expected)


POSITIONAL_COMMAND = "#if $test:\n" "    $test\n" "#end if\n"

NON_POSITIONAL_NON_FLAG_TEXT = "#if $test:\n" "    --test '$test'\n" "#end if\n"

FLAG_COMMAND = "$test\n"

FLAG_COMMAND_TEXT = "$test\n"

REPEAT_COMMAND_TEXT_DATA = (
    "#for $item in $test:\n" "    #if $item:\n" "        --test '$item'\n" "    #end if\n" "#end for\n"
)

REPEAT_FLAG_COMMAND = "#for $item in $test:\n" "    $item\n" "#end for\n"


def test_no_action_positional_command():
    command = extract_command("none_positional")
    assert_equal(command, POSITIONAL_COMMAND)


def test_no_action_command():
    command = extract_command("none")
    assert_equal(command, NON_POSITIONAL_NON_FLAG_TEXT)


def test_store_command():
    command = extract_command("store")
    assert_equal(command, NON_POSITIONAL_NON_FLAG_TEXT)


def test_store_const_command():
    command = extract_command("store_const")
    assert_equal(command, FLAG_COMMAND)


def test_store_const_command_text_data():
    command = extract_command("store_const_text")
    assert_equal(command, FLAG_COMMAND)


def test_store_true_command():
    command = extract_command("store_true")
    assert_equal(command, FLAG_COMMAND)


def test_append_command():
    command = extract_command("append")
    assert_equal(command, REPEAT_COMMAND_TEXT_DATA)


def test_append_const_command():
    command = extract_command("append_const")

    assert_equal(command, REPEAT_FLAG_COMMAND)


def test_count_command():
    command = extract_command("count")

    assert_equal(command, REPEAT_FLAG_COMMAND)


def test_version_command():
    command = extract_command("version")
    assert_equal(command, "")


@pytest.mark.skipif(sys.version_info < (3, 8), reason="Extend command is supported from python 3.8")
def test_extend_command():
    command = extract_command("extend")
    assert_equal(command, REPEAT_COMMAND_TEXT_DATA)


def extract_xml(func_name: str) -> Tuple[str, str, Optional[ParamInfo]]:
    module = _prepare_test_module(func_name)
    parser = obtain_parser(module)

    if parser is None:
        return "", "", None

    data_inputs: Dict[str, str] = dict()
    reserved_names: Set[str] = set()
    name_map: Dict[str, str] = dict()
    section_map: Dict[str, str] = dict()

    auto_inputs, auto_outputs, version_command_param = xml_from_decoy(
        parser, data_inputs, reserved_names, name_map, section_map
    )

    return xml_to_string(auto_inputs, 0), xml_to_string(auto_outputs, 0), version_command_param


def extract_command(func_name: str) -> str:
    module = _prepare_test_module(func_name)
    parser = obtain_parser(module)

    data_inputs: Dict[str, str] = dict()
    reserved_names: Set[str] = set()
    name_map: Dict[str, str] = dict()
    section_map: Dict[str, str] = dict()

    if parser is None:
        return ""

    return command_from_decoy(parser, data_inputs, reserved_names, name_map, section_map, skip_default_namespace=True)


def _prepare_test_module(func_name: str) -> ast.Module:
    return load_function_body(os.path.join(TEST_AUTOPYGEN_DATA, "autopygen_actions_cases.py"), func_name)
