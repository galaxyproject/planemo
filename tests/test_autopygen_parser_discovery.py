import ast
import os

from planemo.autopygen.source_file_parsing.decoy_parser import obtain_class_def
from planemo.autopygen.source_file_parsing.parser_discovery_and_init import (
    ArgumentCreationDiscovery,
    GroupAndSubparsersDiscovery,
    ImportDiscovery,
    SimpleParserDiscoveryAndReplacement,
)
from .test_utils import (
    assert_equal,
    load_function_body,
    TEST_AUTOPYGEN_DATA,
)


def test_import_discovery():
    result = []
    expected_known_names = {"asd", "parser"}
    discovery = ImportDiscovery(result)
    module = _prepare_test_module("import_discovery")

    actions, argparse_module_alias, argparse_class_alias, known_names = discovery.visit_and_report(module)

    assert_equal(len(actions), 4)
    assert_equal(argparse_module_alias, "parser")
    assert_equal(known_names, expected_known_names)


def test_parser_discovery_and_replacement():
    result = []
    module = _prepare_test_module("parser_discovery_and_replacement")
    discovery = ImportDiscovery(result)
    custom_parser_class_def = obtain_class_def()

    actions, argparse_module_alias, argparse_class_alias, known_names = discovery.visit_and_report(module)

    discovery = SimpleParserDiscoveryAndReplacement(
        actions, argparse_module_alias, argparse_class_alias, custom_parser_class_def
    )

    actions, parser_name = discovery.visit_and_report(module)

    assert_equal(len(actions), 3)
    _, decoy_def, parser_init = actions
    assert_equal(type(decoy_def), ast.ClassDef)
    assert_equal(decoy_def.name, "DecoyParser")
    assert_equal(type(parser_init), ast.Assign)
    assert_equal(parser_init.value.func.id, "DecoyParser")


def test_group_discovery():
    result = []

    module = _prepare_test_module("group_discovery")

    discovery = ImportDiscovery(result)
    custom_parser_class_def = obtain_class_def()

    actions, argparse_module_alias, argparse_class_alias, known_names = discovery.visit_and_report(module)

    discovery = SimpleParserDiscoveryAndReplacement(
        actions, argparse_module_alias, argparse_class_alias, custom_parser_class_def
    )

    actions, parser_name = discovery.visit_and_report(module)
    groups = GroupAndSubparsersDiscovery(actions, known_names, parser_name)
    actions, known_names = groups.visit_and_report(module)

    assert_equal(len(actions), 4)
    assert_equal(type(actions[3]), ast.Assign)
    group = actions[3]
    assert_equal(group.targets[0].id, "group")
    assert_equal(group.value.func.attr, "add_argument_group")


def test_argument_creation_discovery():
    result = []
    module = _prepare_test_module("argument_creation")

    discovery = ImportDiscovery(result)
    custom_parser_class_def = obtain_class_def()

    actions, argparse_module_alias, argparse_class_alias, known_names = discovery.visit_and_report(module)

    discovery = SimpleParserDiscoveryAndReplacement(
        actions, argparse_module_alias, argparse_class_alias, custom_parser_class_def
    )

    actions, parser_name = discovery.visit_and_report(module)
    groups = GroupAndSubparsersDiscovery(actions, known_names, parser_name)
    actions, known_names = groups.visit_and_report(module)

    argument_creation = ArgumentCreationDiscovery(actions, parser_name)
    (actions,) = argument_creation.visit_and_report(module)

    assert_equal(len(actions), 4)
    assert_equal(type(actions[3]), ast.Expr)
    arg_creat_call = actions[3]
    assert_equal(arg_creat_call.value.func.attr, "add_argument")

    arg = arg_creat_call.value.args[0]
    if isinstance(arg, ast.Str):
        assert_equal(arg.s, "--test")
    # different python versions produce different results
    else:
        assert_equal(arg.value, "--test")


def _prepare_test_module(func_name: str) -> ast.Module:
    return load_function_body(os.path.join(TEST_AUTOPYGEN_DATA, "autopygen_parser_extraction_test_data.py"), func_name)
