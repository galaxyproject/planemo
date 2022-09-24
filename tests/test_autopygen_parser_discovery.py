import ast
import os

from planemo.autopygen.source_file_parsing.parser_discovery_and_init import SimpleParserDiscoveryAndReplacement, \
    ImportDiscovery, GroupAndSubparsersDiscovery, ArgumentCreationDiscovery, get_parser_init_and_actions
from .test_utils import assert_equal, TEST_DATA_DIR


def test_import_discovery():
    result = []
    expected_actions = []
    expected_names = []

    discovery = ImportDiscovery(result)
    module = _prepare_test_module("test_data_visitors")

    actions, known_names = discovery.visit_and_report(module)

    assert_equal(actions, expected_actions)
    assert_equal(known_names, expected_names)


def test_parser_discovery_and_replacement():
    result = []
    expected_actions = []
    expected_names = []

    discovery = SimpleParserDiscoveryAndReplacement(result)
    module = _prepare_test_module("test_data_visitors")

    actions, known_names = discovery.visit_and_report(module)

    assert_equal(actions, expected_actions)
    assert_equal(known_names, expected_names)


def test_group_discovery():
    result = []
    expected_actions = []
    expected_names = []

    discovery = GroupAndSubparsersDiscovery(result)
    module = _prepare_test_module("test_data_visitors")

    actions, known_names = discovery.visit_and_report(module)

    assert_equal(actions, expected_actions)
    assert_equal(known_names, expected_names)


def test_argument_creation_discovery():
    result = []
    expected_actions = []
    expected_names = []

    discovery = ArgumentCreationDiscovery(result)
    module = _prepare_test_module("test_data_visitors")

    actions, known_names = discovery.visit_and_report(module)

    assert_equal(actions, expected_actions)
    assert_equal(known_names, expected_names)


def test_simple_file():
    """Test very simple command."""
    SimpleParserDiscoveryAndReplacement()
    assert_equal(command_io.cheetah_template, "random_fastq")


def test_parser_in_function():
    """Test example input/output transition."""
    command_io = _example("convert 1.bed 1.bam", example_outputs=["1.bam"], example_inputs=["1.bed"])
    assert_equal(command_io.cheetah_template, "convert '$input1' '$output1'")
    assert_equal(len(command_io.outputs), 1)
    assert_equal(len(command_io.inputs), 1)


def test_unknown_import():
    """Test example input/output transition."""
    command_io = _example('convert "1.bed" "1.bam"', example_outputs=["1.bam"], example_inputs=["1.bed"])
    assert_equal(command_io.cheetah_template, "convert '$input1' '$output1'")


def test_no_parser_found():
    """Test example input/output transition."""
    command_io = _example("convert '1.bed' '1.bam'", example_outputs=["1.bam"], example_inputs=["1.bed"])
    assert_equal(command_io.cheetah_template, "convert '$input1' '$output1'")


def test_custom_parser_used():
    command_io = _example("convert '1.bed' '1.bam'", example_outputs=["1.bam"], example_inputs=["1.bed"])
    cwl_properties = command_io.cwl_properties()

    assert_equal(cwl_properties["base_command"], ["convert"])
    assert_equal(cwl_properties["inputs"][0].position, 1)
    assert_equal(cwl_properties["outputs"][0].position, 2)


def test_subparsers():
    command_io = _example("convert '1.bed' '1.bam'", example_outputs=["1.bam"], example_inputs=["1.bed"])
    cwl_properties = command_io.cwl_properties()

    assert_equal(cwl_properties["base_command"], ["convert"])
    assert_equal(cwl_properties["inputs"][0].position, 1)
    assert_equal(cwl_properties["outputs"][0].position, 2)


def test_argument_groups():
    command_io = _example("convert '1.bed' '1.bam'", example_outputs=["1.bam"], example_inputs=["1.bed"])
    cwl_properties = command_io.cwl_properties()

    assert_equal(cwl_properties["base_command"], ["convert"])
    assert_equal(cwl_properties["inputs"][0].position, 1)
    assert_equal(cwl_properties["outputs"][0].position, 2)


def _prepare_test_module(func_name: str) -> ast.Module:
    test_file = open(os.join(TEST_DATA_DIR, ))
