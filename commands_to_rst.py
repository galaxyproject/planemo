#!/usr/bin/env python

import os
from string import Template

from click.testing import CliRunner
from planemo.cli import list_cmds
from planemo import cli

planemo_cli = cli.planemo
runner = CliRunner()

INTERNAL_COMMANDS = ['travis_before_install']

COMMAND_TEMPLATE = Template('''
``${command}`` command
===============================

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
        if line.startswith("    "):
            return line[4:]
        else:
            return line
    clean_rst = "\n".join(map(clean_rst_line, raw_rst.split("\n")))

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
        elif option_lines:
            new_lines.append("    %s" % line)
    text = COMMAND_TEMPLATE.safe_substitute(
        command=command,
        command_help="\n".join(new_lines),
    )
    commands += "\n.. include:: commands/%s.rst" % command
    open(os.path.join(command_doc_dir, command + ".rst"), "w").write(text)

open(os.path.join("docs", "commands.rst"), "w").write(commands)
