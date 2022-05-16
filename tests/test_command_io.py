"""Test cases for the CommandIO abstraction in tool_builder."""

from planemo.tool_builder import CommandIO
from .test_utils import assert_equal


def test_simplest_command():
    """Test very simple command."""
    command_io = _example("random_fastq")
    assert_equal(command_io.cheetah_template, "random_fastq")


def test_example_and_quotes():
    """Test example input/output transition."""
    command_io = _example("convert 1.bed 1.bam", example_outputs=["1.bam"], example_inputs=["1.bed"])
    assert_equal(command_io.cheetah_template, "convert '$input1' '$output1'")
    assert_equal(len(command_io.outputs), 1)
    assert_equal(len(command_io.inputs), 1)


def test_example_already_quoted():
    """Test example input/output transition."""
    command_io = _example('convert "1.bed" "1.bam"', example_outputs=["1.bam"], example_inputs=["1.bed"])
    assert_equal(command_io.cheetah_template, "convert '$input1' '$output1'")


def test_example_already_quoted_single():
    """Test example input/output transition."""
    command_io = _example("convert '1.bed' '1.bam'", example_outputs=["1.bam"], example_inputs=["1.bed"])
    assert_equal(command_io.cheetah_template, "convert '$input1' '$output1'")


def test_example_cwl_simplest():
    command_io = _example("convert '1.bed' '1.bam'", example_outputs=["1.bam"], example_inputs=["1.bed"])
    cwl_properties = command_io.cwl_properties()

    assert_equal(cwl_properties["base_command"], ["convert"])
    assert_equal(cwl_properties["inputs"][0].position, 1)
    assert_equal(cwl_properties["outputs"][0].position, 2)


def test_example_cwl_argument():
    command_io = _example("seqtk convert '1.bed' moo '1.bam'", example_outputs=["1.bam"], example_inputs=["1.bed"])
    cwl_properties = command_io.cwl_properties()

    assert_equal(cwl_properties["base_command"], ["seqtk", "convert"])
    assert_equal(cwl_properties["inputs"][0].position, 1)
    assert_equal(cwl_properties["inputs"][0].prefix, None)
    assert_equal(cwl_properties["outputs"][0].position, 3)
    assert_equal(cwl_properties["outputs"][0].prefix, None)
    assert_equal(cwl_properties["outputs"][0].glob, "$(inputs.output1)")
    assert_equal(len(cwl_properties["arguments"]), 1)
    assert_equal(cwl_properties["arguments"][0].position, 2)
    assert_equal(cwl_properties["arguments"][0].value, "moo")


def test_example_cwl_simple_redirect():
    command_io = _example("seqtk convert '1.bed' > '1.bam'", example_outputs=["1.bam"], example_inputs=["1.bed"])
    cwl_properties = command_io.cwl_properties()

    assert_equal(cwl_properties["base_command"], ["seqtk", "convert"])
    assert_equal(cwl_properties["inputs"][0].position, 1)
    assert_equal(cwl_properties["outputs"][0].glob, "out")
    assert_equal(cwl_properties["stdout"], "out")


def test_prefixes_separated():
    command_io = _example(
        "seqtk convert -i '1.bed' --output '1.bam'", example_outputs=["1.bam"], example_inputs=["1.bed"]
    )
    cwl_properties = command_io.cwl_properties()
    assert_equal(cwl_properties["base_command"], ["seqtk", "convert"])
    assert_equal(cwl_properties["inputs"][0].position, 1)
    assert_equal(cwl_properties["inputs"][0].prefix.prefix, "-i")
    assert_equal(cwl_properties["inputs"][0].prefix.separated, True)

    assert_equal(cwl_properties["outputs"][0].glob, "$(inputs.output1)")
    assert_equal(cwl_properties["outputs"][0].prefix.prefix, "--output")
    assert_equal(cwl_properties["outputs"][0].prefix.separated, True)
    assert_equal(cwl_properties["stdout"], None)


def test_prefixes_joined():
    command_io = _example("seqtk convert INPUT=1.bed OUTPUT=1.bam", example_outputs=["1.bam"], example_inputs=["1.bed"])
    cwl_properties = command_io.cwl_properties()
    assert_equal(cwl_properties["base_command"], ["seqtk", "convert"])
    assert_equal(cwl_properties["inputs"][0].position, 1)
    assert_equal(cwl_properties["inputs"][0].prefix.prefix, "INPUT=")
    assert_equal(cwl_properties["inputs"][0].prefix.separated, False)

    assert_equal(cwl_properties["outputs"][0].glob, "$(inputs.output1)")
    assert_equal(cwl_properties["outputs"][0].prefix.prefix, "OUTPUT=")
    assert_equal(cwl_properties["outputs"][0].prefix.separated, False)
    assert_equal(cwl_properties["stdout"], None)


def test_integer_parameters():
    command_io = _example(
        "seqtk convert --size 100 -i '1.bed' --threshold 2.0 --output_type bam > '1.bam'",
        example_outputs=["1.bam"],
        example_inputs=["1.bed"],
    )
    cwl_properties = command_io.cwl_properties()
    assert_equal(cwl_properties["base_command"], ["seqtk", "convert"])
    assert_equal(len(cwl_properties["inputs"]), 4)
    assert_equal(cwl_properties["inputs"][0].position, 1)
    assert_equal(cwl_properties["inputs"][0].type, "int")
    assert_equal(cwl_properties["inputs"][0].prefix.prefix, "--size")

    assert_equal(cwl_properties["inputs"][1].position, 2)
    assert_equal(cwl_properties["inputs"][1].type, "File")
    assert_equal(cwl_properties["inputs"][1].prefix.prefix, "-i")

    assert_equal(cwl_properties["inputs"][2].position, 3)
    assert_equal(cwl_properties["inputs"][2].type, "float")
    assert_equal(cwl_properties["inputs"][2].prefix.prefix, "--threshold")

    assert_equal(cwl_properties["inputs"][3].position, 4)
    assert_equal(cwl_properties["inputs"][3].type, "string")
    assert_equal(cwl_properties["inputs"][3].prefix.prefix, "--output_type")

    assert_equal(cwl_properties["outputs"][0].glob, "out")
    assert_equal(cwl_properties["stdout"], "out")


def _example(example_command, example_outputs=[], example_inputs=[]):
    """Build a CommandIO object for test cases."""
    kwds = {}
    kwds["example_command"] = example_command
    kwds["example_output"] = example_outputs
    kwds["example_input"] = example_inputs
    return CommandIO(**kwds)
