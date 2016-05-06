"""Test cases for the CommandIO abstraction in tool_builder."""

from planemo.tool_builder import CommandIO


def test_simplest_command():
    """Test very simple command."""
    command_io = _example("random_fastq")
    _assert_eq(command_io.cheetah_template, "random_fastq")


def test_example_and_quotes():
    """Test example input/output transition."""
    command_io = _example("convert 1.bed 1.bam", example_outputs=["1.bam"], example_inputs=["1.bed"])
    _assert_eq(command_io.cheetah_template, 'convert "$input1" "$output1"')
    _assert_eq(len(command_io.outputs), 1)
    _assert_eq(len(command_io.inputs), 1)


def test_example_already_quoted():
    """Test example input/output transition."""
    command_io = _example('convert "1.bed" "1.bam"', example_outputs=["1.bam"], example_inputs=["1.bed"])
    _assert_eq(command_io.cheetah_template, 'convert "$input1" "$output1"')


def test_example_already_quoted_single():
    """Test example input/output transition."""
    command_io = _example("convert '1.bed' '1.bam'", example_outputs=["1.bam"], example_inputs=["1.bed"])
    _assert_eq(command_io.cheetah_template, "convert '$input1' '$output1'")


def test_example_cwl_simplest():
    command_io = _example("convert '1.bed' '1.bam'", example_outputs=["1.bam"], example_inputs=["1.bed"])
    cwl_properties = command_io.cwl_properties()

    _assert_eq(cwl_properties["base_command"], ["convert"])
    _assert_eq(cwl_properties["inputs"][0].position, 1)
    _assert_eq(cwl_properties["outputs"][0].position, 2)


def test_example_cwl_argument():
    command_io = _example("seqtk convert '1.bed' moo '1.bam'", example_outputs=["1.bam"], example_inputs=["1.bed"])
    cwl_properties = command_io.cwl_properties()

    _assert_eq(cwl_properties["base_command"], ["seqtk", "convert"])
    _assert_eq(cwl_properties["inputs"][0].position, 1)
    _assert_eq(cwl_properties["inputs"][0].prefix, None)
    _assert_eq(cwl_properties["outputs"][0].position, 3)
    _assert_eq(cwl_properties["outputs"][0].prefix, None)
    _assert_eq(cwl_properties["outputs"][0].glob, "$(inputs.output1)")
    _assert_eq(len(cwl_properties["arguments"]), 1)
    _assert_eq(cwl_properties["arguments"][0].position, 2)
    _assert_eq(cwl_properties["arguments"][0].value, "moo")


def test_example_cwl_simple_redirect():
    command_io = _example("seqtk convert '1.bed' > '1.bam'", example_outputs=["1.bam"], example_inputs=["1.bed"])
    cwl_properties = command_io.cwl_properties()

    _assert_eq(cwl_properties["base_command"], ["seqtk", "convert"])
    _assert_eq(cwl_properties["inputs"][0].position, 1)
    _assert_eq(cwl_properties["outputs"][0].glob, "out")
    _assert_eq(cwl_properties["stdout"], "out")


def test_prefixes_separated():
    command_io = _example("seqtk convert -i '1.bed' --output '1.bam'", example_outputs=["1.bam"], example_inputs=["1.bed"])
    cwl_properties = command_io.cwl_properties()
    _assert_eq(cwl_properties["base_command"], ["seqtk", "convert"])
    _assert_eq(cwl_properties["inputs"][0].position, 1)
    _assert_eq(cwl_properties["inputs"][0].prefix.prefix, "-i")
    _assert_eq(cwl_properties["inputs"][0].prefix.separated, True)

    _assert_eq(cwl_properties["outputs"][0].glob, "$(inputs.output1)")
    _assert_eq(cwl_properties["outputs"][0].prefix.prefix, "--output")
    _assert_eq(cwl_properties["outputs"][0].prefix.separated, True)
    _assert_eq(cwl_properties["stdout"], None)


def test_prefixes_joined():
    command_io = _example("seqtk convert INPUT=1.bed OUTPUT=1.bam", example_outputs=["1.bam"], example_inputs=["1.bed"])
    cwl_properties = command_io.cwl_properties()
    _assert_eq(cwl_properties["base_command"], ["seqtk", "convert"])
    _assert_eq(cwl_properties["inputs"][0].position, 1)
    _assert_eq(cwl_properties["inputs"][0].prefix.prefix, "INPUT=")
    _assert_eq(cwl_properties["inputs"][0].prefix.separated, False)

    _assert_eq(cwl_properties["outputs"][0].glob, "$(inputs.output1)")
    _assert_eq(cwl_properties["outputs"][0].prefix.prefix, "OUTPUT=")
    _assert_eq(cwl_properties["outputs"][0].prefix.separated, False)
    _assert_eq(cwl_properties["stdout"], None)


def _example(example_command, example_outputs=[], example_inputs=[]):
    """Build a CommandIO object for test cases."""
    kwds = {}
    kwds["example_command"] = example_command
    kwds["example_output"] = example_outputs
    kwds["example_input"] = example_inputs
    return CommandIO(**kwds)


def _assert_eq(a, b):
    assert a == b, "%s != %s" % (a, b)
