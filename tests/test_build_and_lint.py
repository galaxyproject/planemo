
from .test_utils import CliTestCase


INIT_COMMAND = [
    "tool_init", "--force",
    "--id", "seqtk_seq",
    "--name", "Convert to FASTA (seqtk)",
    "--requirement", "seqtk@1.0-r68",
    "--example_command", "seqtk seq -a 2.fastq > 2.fasta",
    "--example_input", "2.fastq",
    "--example_output", "2.fasta",
    "--test_case",
    "--help_text", "The help text.",
    "--doi", "10.1101/014043"
]


class BuildAndLintTestCase(CliTestCase):

    def test_build_and_lint(self):
        with self._isolate():
            self._check_exit_code(_init_command())
            self._check_lint(exit_code=0)

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

    def _check_lint(self, exit_code=0):
        lint_cmd = ["lint", "--fail_level", "warn", "seqtk_seq.xml"]
        self._check_exit_code(lint_cmd, exit_code=exit_code)


def _init_command(test_case=True, help_text=True, doi=True):
    command = [
        "tool_init", "--force",
        "--id", "seqtk_seq",
        "--name", "Convert to FASTA (seqtk)",
        "--requirement", "seqtk@1.0-r68",
        "--example_command", "seqtk seq -a 2.fastq > 2.fasta",
        "--example_input", "2.fastq",
        "--example_output", "2.fasta"
    ]
    if test_case:
        command.append("--test_case")
    if help_text:
        command.extend(["--help_text", "The help text."])
    if doi:
        command.extend(["--doi", "10.1101/014043"])
    return command
