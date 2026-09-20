import os

import yaml

from .test_utils import (
    CliTestCase,
    skip_if_environ,
)


class BuildAndLintTestCase(CliTestCase):
    def test_build_and_lint(self):
        with self._isolate():
            self._check_exit_code(_init_command())
            self._check_lint(exit_code=0)

    def test_build_and_lint_with_macros(self):
        with self._isolate() as f:
            self._check_exit_code(_init_command(macros=True))
            self._check_lint(exit_code=0)
            macros_file = os.path.join(f, "macros.xml")
            assert os.path.exists(macros_file)

    def test_lint_fails_if_no_help(self):
        with self._isolate():
            self._check_exit_code(_init_command(help_text=False))
            self._check_lint(exit_code=1)

    def test_lint_fails_if_no_test(self):
        with self._isolate():
            self._check_exit_code(_init_command(test_case=False))
            self._check_lint(exit_code=1)

    def test_lint_fails_if_no_doi(self):
        with self._isolate():
            self._check_exit_code(_init_command(doi=False))
            self._check_lint(exit_code=1)

    @skip_if_environ("PLANEMO_SKIP_CWLTOOL_TESTS")
    def test_cwl(self):
        with self._isolate() as f:
            self._check_exit_code(_cwl_init_command())
            self._check_lint(filename="seqtk_seq.cwl", exit_code=0)

            with open(os.path.join(f, "seqtk_seq.cwl")) as stream:
                process_dict = yaml.safe_load(stream)
            assert process_dict["id"] == "seqtk_seq"
            assert process_dict["label"] == "Convert to FASTA (seqtk)"
            assert process_dict["baseCommand"] == ["seqtk", "seq"]
            input0 = process_dict["inputs"]["input1"]
            assert input0["inputBinding"]["position"] == 1
            assert input0["inputBinding"]["prefix"] == "-A"
            assert input0["type"] == "File"
            output = process_dict["outputs"]["output1"]
            assert output["type"] == "File"
            assert output["outputBinding"]["glob"] == "out"
            assert process_dict["stdout"] == "out"

            with open(os.path.join(f, "seqtk_seq_tests.yml")) as stream:
                test_dict = yaml.safe_load(stream)
                assert test_dict

    @skip_if_environ("PLANEMO_SKIP_CWLTOOL_TESTS")
    def test_cwl_fail_on_empty_help(self):
        with self._isolate():
            self._check_exit_code(_cwl_init_command(help_text=False))
            self._check_lint(filename="seqtk_seq.cwl", exit_code=1)

    @skip_if_environ("PLANEMO_SKIP_CWLTOOL_TESTS")
    def test_cwl_fail_on_no_docker(self):
        with self._isolate():
            self._check_exit_code(_cwl_init_command(help_text=False))
            self._check_lint(filename="seqtk_seq.cwl", exit_code=1)

    def _check_lint(self, filename="seqtk_seq.xml", exit_code=0):
        lint_cmd = ["lint", "--fail_level", "warn", filename]
        try:
            self._check_exit_code(lint_cmd, exit_code=exit_code)
        except Exception:
            with open(filename) as f:
                print("Failing file contents are [%s]." % f.read())
            raise


def _cwl_init_command(help_text=True, container=True, test_case=True):
    command = [
        "tool_init",
        "--force",
        "--cwl",
        "--id",
        "seqtk_seq",
        "--name",
        "Convert to FASTA (seqtk)",
        "--name",
        "Convert to FASTA (seqtk)",
        "--example_command",
        "seqtk seq -A 2.fastq > 2.fasta",
        "--example_input",
        "2.fastq",
        "--example_output",
        "2.fasta",
    ]
    if container:
        command.extend(["--container", "quay.io/biocontainers/seqtk:1.2--0"])
    if help_text:
        command.extend(["--help_text", "The help text."])
    if test_case:
        command.append("--test_case")

    return command


def _init_command(test_case=True, help_text=True, doi=True, macros=False):
    command = [
        "tool_init",
        "--force",
        "--id",
        "seqtk_seq",
        "--name",
        "Convert to FASTA (seqtk)",
        "--requirement",
        "seqtk@1.0-r68",
        "--example_command",
        "seqtk seq -A 2.fastq > 2.fasta",
        "--example_input",
        "2.fastq",
        "--example_output",
        "2.fasta",
    ]
    if test_case:
        command.append("--test_case")
    if help_text:
        command.extend(["--help_text", "The help text."])
    if doi:
        command.extend(["--doi", "10.1101/014043"])
        command.extend(["--cite_url", "https://github.com/ekg/vcflib"])
        command.extend(["--cite_url", "http://wiki.hpc.ufl.edu/doc/Seqtk"])
    if macros:
        command.append("--macros")
    return command
