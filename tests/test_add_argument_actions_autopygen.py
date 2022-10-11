import ast
import os
from typing import Tuple, Optional

from planemo.autopygen.argument_parser_conversion import command_from_decoy, xml_from_decoy, \
    obtain_parser
from planemo.autopygen.param_info import ParamInfo
from tests.test_utils import load_function_body, TEST_AUTOPYGEN_DATA, assert_equal


def test_no_action_positional():
    inputs, _, _ = extract_xml("none_positional")
    assert_equal(inputs, '<param argument="test" type="text" optional="false" label="test" />\n')


def test_no_action():
    inputs, _, _ = extract_xml("none")
    assert_equal(inputs, '<param argument="--test" type="text" optional="true" label="test" />\n')


def test_store():
    inputs, _, _ = extract_xml("store")
    assert_equal(inputs, '<param argument="--test" type="text" optional="true" label="test" />\n')


def test_store_const():
    inputs, _, _ = extract_xml("store_const")
    assert_equal(inputs, '<param argument="--test" type="boolean" truevalue="--test" falsevalue="" checked="false" '
                         'optional="true" label="test" />\n')


def test_store_true():
    inputs, _, _ = extract_xml("store_true")
    assert_equal(inputs, '<param argument="--test" type="boolean" truevalue="--test" falsevalue="" checked="false" '
                         'optional="true" label="test" />\n')


def test_append():
    inputs, _, _ = extract_xml("append")
    expected = ('<repeat name="test_repeat" title="test_repeat">\n'
                '    <param argument="--test" type="text" optional="true" label="test" />\n'
                '</repeat>\n')

    assert_equal(inputs, expected)


def test_append_const():
    inputs, _, _ = extract_xml("append_const")
    expected = ('<repeat name="test_repeat" title="test_repeat">\n'
                '    <param argument="--test" type="boolean" truevalue="--test" falsevalue="" checked="false" '
                'optional="true" label="test" />\n'
                '</repeat>\n')
    assert_equal(inputs, expected)


def test_count():
    inputs, _, _ = extract_xml("count")
    expected = ('<repeat name="test_repeat" title="test_repeat">\n'
                '    <param argument="--test" type="boolean" truevalue="--test" falsevalue="" checked="false" '
                'optional="true" label="test" />\n'
                '</repeat>\n')
    assert_equal(inputs, expected)


def test_version():
    inputs, _, param_info = extract_xml("version")

    assert_equal(inputs, '')
    assert_equal(type(param_info), ParamInfo)
    assert_equal(param_info.param_type.is_version, True)


def test_extend():
    inputs, _, _ = extract_xml("extend")
    expected = ('<repeat name="test_repeat" title="test_repeat">\n'
                '    <param argument="--test" type="text" optional="true" label="test" />\n'
                '</repeat>\n')
    assert_equal(inputs, expected)


POSITIONAL_COMMAND = (
    '## test definition\n'
    '#if $test:\n'
    '   $test\n'
    '#end if\n'
    '## end test definition\n'
)

NON_POSITIONAL_NON_FLAG = (
    '## test definition\n'
    '#if $test:\n'
    '   --test $test\n'
    '#end if\n'
    '## end test definition\n'
)

FLAG_COMMAND = (
    '## FLAG test definition\n'
    '$test\n'
    '## end FLAG test definition\n'
)

REPEAT_COMMAND = (
    '## test definition\n'
    '#for $item in $test:\n'
    '   ## test definition\n'
    '   #if $item:\n'
    '      --test $item\n'
    '   #end if\n'
    '   ## end test definition\n'
    '\n'
    '#end for\n'
    '## end test definition\n'
)

REPEAT_FLAG_COMMAND = (
    '## test definition\n'
    '#for $item in $test:\n'
    '   ## FLAG test definition\n'
    '   $item\n'
    '   ## end FLAG test definition\n'
    '\n'
    '#end for\n'
    '## end test definition\n'
)


def test_no_action_positional_command():
    command = extract_command("none_positional")
    assert_equal(command, POSITIONAL_COMMAND)


def test_no_action_command():
    command = extract_command("none")
    assert_equal(command, NON_POSITIONAL_NON_FLAG)


def test_store_command():
    command = extract_command("store")
    assert_equal(command, NON_POSITIONAL_NON_FLAG)


def test_store_const_command():
    command = extract_command("store_const")
    assert_equal(command, FLAG_COMMAND)


def test_store_true_command():
    command = extract_command("store_true")
    assert_equal(command, FLAG_COMMAND)


def test_append_command():
    command = extract_command("append")
    assert_equal(command, REPEAT_COMMAND)


def test_append_const_command():
    command = extract_command("append_const")

    assert_equal(command, REPEAT_FLAG_COMMAND)


def test_count_command():
    command = extract_command("count")

    assert_equal(command, REPEAT_FLAG_COMMAND)


def test_version_command():
    command = extract_command("version")
    assert_equal(command, '')


def test_extend_command():
    command = extract_command("extend")
    assert_equal(command, REPEAT_COMMAND)


def extract_xml(func_name: str) -> Tuple[str, str, Optional[ParamInfo]]:
    module = _prepare_test_module(func_name)
    parser = obtain_parser(module)

    data_inputs = dict()
    reserved_names = set()
    name_map = dict()
    section_map = dict()
    auto_inputs, auto_outputs, version_command_param = \
        xml_from_decoy(parser, data_inputs, reserved_names, name_map, section_map, depth=0)

    return auto_inputs, auto_outputs, version_command_param


def extract_command(func_name: str) -> str:
    module = _prepare_test_module(func_name)
    parser = obtain_parser(module)

    data_inputs = dict()
    reserved_names = set()
    name_map = dict()
    section_map = dict()

    return command_from_decoy(parser, data_inputs, reserved_names, name_map, section_map, skip_default_namespace=True)


def _prepare_test_module(func_name: str) -> ast.Module:
    return load_function_body(os.path.join(TEST_AUTOPYGEN_DATA, "autopygen_actions_cases.py"), func_name)
