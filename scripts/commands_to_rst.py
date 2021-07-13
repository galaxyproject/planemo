#!/usr/bin/env python

import os
import sys
from string import Template

from click.testing import CliRunner

project_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, project_dir)

from planemo import cli  # noqa: E402
from planemo.cli import list_cmds  # noqa: E402

planemo_cli = cli.planemo
runner = CliRunner()

# Don't document the following commands - they should not be considered part
# of the planemo API.
INTERNAL_COMMANDS = [
    'create_gist',
    'shed_download',
]

COMMAND_TEMPLATE = Template('''
``${command}`` command
======================================

This section is auto-generated from the help text for the planemo command
``${command}``. This help message can be generated with ``planemo ${command}
--help``.

${command_help}
''')

COMMANDS_TEMPLATE = """========
Commands
========

Planemo is a set of utilities for developing Galaxy tools. Each utility is
implemented as a subcommand of the ``planemo`` executable. This section of the
documentation describes these commands.

"""

command_doc_dir = os.path.join("docs", "commands")
commands = COMMANDS_TEMPLATE

for command in list_cmds():
    if command in INTERNAL_COMMANDS:
        continue

    command_obj = cli.name_to_command(command)
    function = command_obj.callback
    raw_rst = function.__doc__

    def clean_rst_line(line):
        # Check for that click \b.
        if line.startswith('    \x08'):
            return "::\n"
        if line.startswith("    "):
            return line[4:]
        else:
            return line
    all_lines = raw_rst.split("\n")
    all_lines[0] = all_lines[0].lstrip()  # """ Fix docs like this.
    clean_rst = "\n".join(map(clean_rst_line, all_lines))

    result = runner.invoke(planemo_cli, [command, "--help"])
    output = result.output
    lines = output.split("\n")
    new_lines = []
    help_lines = False
    option_lines = False
    for line in lines:
        if line.startswith("Usage: "):
            new_lines.append("**Usage**::\n\n    %s" % line[len("Usage: "):])
            new_lines.append("\n**Help**\n")
            new_lines.append(clean_rst)
            help_lines = True
        elif line.startswith("Options:"):
            help_lines = False
            new_lines.append("**Options**::\n\n")
            option_lines = True
        elif line.startswith('\b') or line.startswith("^H") or line.startswith('\x08') or line.startswith(''):
            new_lines.append("\n::\n\n")
        elif option_lines:
            new_lines.append("    %s" % line)
    text = COMMAND_TEMPLATE.safe_substitute(
        command=command,
        command_help="\n".join(new_lines),
    )
    commands += "\n.. include:: commands/%s.rst" % command
    with open(os.path.join(command_doc_dir, command + ".rst"), "w") as fh:
        fh.write(text)

with open(os.path.join("docs", "commands.rst"), "w") as fh:
    fh.write(commands)
