"""This module contains :func:`build` to build tool descriptions.

This class is used by the `tool_init` command and can be used to build
Galaxy and CWL tool descriptions.
"""
from collections import namedtuple
import re
import shlex
import subprocess

from planemo import templates
from planemo import io


TOOL_TEMPLATE = """<tool id="{{id}}" name="{{name}}" version="{{version}}">
{%- if description %}
    <description>{{ description }}</description>
{%- endif %}
{%- if macros %}
    <macros>
        <import>macros.xml</import>
    </macros>
    <expand macro="requirements" />
    <expand macro="stdio" />
{%- if version_command %}
    <expand macro="version_command" />
{%- endif %}
{%- else %}
    <requirements>
{%- for requirement in requirements %}
        {{ requirement }}
{%- endfor %}
{%- for container in containers %}
        {{ container }}
{%- endfor %}
    </requirements>
    <stdio>
        <exit_code range="1:" />
    </stdio>
{%- if version_command %}
    <version_command>{{ version_command }}</version_command>
{%- endif %}
{%- endif %}
    <command><![CDATA[
{%- if command %}
        {{ command }}
{%- else %}
        TODO: Fill in command template.
{%- endif %}
    ]]></command>
    <inputs>
{%- for input in inputs %}
        {{ input }}
{%- endfor %}
    </inputs>
    <outputs>
{%- for output in outputs %}
        {{ output }}
{%- endfor %}
    </outputs>
{%- if tests %}
    <tests>
{%- for test in tests %}
        <test>
{%- for param in test.params %}
            <param name="{{ param[0]}}" value="{{ param[1] }}"/>
{%- endfor %}
{%- for output in test.outputs %}
            <output name="{{ output[0] }}" file="{{ output[1] }}"/>
{%- endfor %}
        </test>
{%- endfor %}
    </tests>
{%- endif %}
    <help><![CDATA[
{%- if help %}
        {{ help }}
{%- else %}
        TODO: Fill in help.
{%- endif %}
    ]]></help>
{%- if macros %}
    <expand macro="citations" />
{%- else %}
{%- if doi or bibtex_citations %}
    <citations>
{%- for single_doi in doi %}
        <citation type="doi">{{ single_doi }}</citation>
{%- endfor %}
{%- for bibtex_citation in bibtex_citations %}
        <citation type="bibtex">{{ bibtex_citation }}</citation>
{%- endfor %}
    </citations>
{%- endif %}
{%- endif %}
</tool>
"""

MACROS_TEMPLATE = """<macros>
    <xml name="requirements">
        <requirements>
{%- for requirement in requirements %}
        {{ requirement }}
{%- endfor %}
            <yield/>
{%- for container in containers %}
        {{ container }}
{%- endfor %}
        </requirements>
    </xml>
    <xml name="stdio">
        <stdio>
            <exit_code range="1:" />
        </stdio>
    </xml>
    <xml name="citations">
        <citations>
{%- for single_doi in doi %}
            <citation type="doi">{{ single_doi }}</citation>
{%- endfor %}
{%- for bibtex_citation in bibtex_citations %}
            <citation type="bibtex">{{ bibtex_citation }}</citation>
{%- endfor %}
            <yield />
        </citations>
    </xml>
{%- if version_command %}
    <xml name="version_command">
        <version_command>{{ version_command }}</version_command>
    </xml>
{%- endif %}
</macros>
"""


CWL_TEMPLATE = """#!/usr/bin/env cwl-runner
cwlVersion: '{{cwl_version}}'
class: CommandLineTool
id: "{{id}}"
label: "{{label}}"
{%- if containers %}
requirements:
{%- for container in containers %}
  - class: DockerRequirement
    dockerPull: {{ container.image_id }}
{%- endfor %}
{%- endif %}
{%- if inputs or outputs %}
inputs:
{%- for input in inputs %}
  - id: {{ input.id }}
    type: {{ input.type }}
    description: |
      TODO
    inputBinding:
      position: {{ input.position }}
{%- if input.prefix %}
      prefix: "{{input.prefix.prefix}}"
{%- if not input.prefix.separated %}
      separate: false
{%- endif %}
{%- endif %}
{%- endfor %}
{%- for output in outputs %}
{%- if output.require_filename %}
  - id: {{ output.id }}
    type: string
    description: |
      Filename for output {{ output.id }}
    inputBinding:
      position: {{ output.position }}
{%- if output.prefix %}
      prefix: "{{output.prefix.prefix}}"
{%- if not output.prefix.separated %}
      separate: false
{%- endif %}
{%- endif %}
{%- endif %}
{%- endfor %}
{%- else %}
inputs: [] # TODO
{%- endif %}
{%- if outputs %}
outputs:
{%- for output in outputs %}
  - id: {{ output.id }}
    type: File
    outputBinding:
      glob: {{ output.glob }}
{%- endfor %}
{%- else %}
outputs: [] # TODO
{%- endif %}
{%- if base_command %}
baseCommand:
{%- for base_command_part in base_command %}
  - "{{ base_command_part}}"
{%- endfor %}
{%- else %}
baseCommand: []
{%- endif %}
{%- if arguments %}
arguments:
{%- for argument in arguments %}
  - valueFrom: "{{ argument.value }}"
    position: {{ argument.position }}
{%- if argument.prefix %}
      prefix: "{{argument.prefix.prefix}}"
{%- if not argument.prefix.separated %}
      separate: false
{%- endif %}
{%- endif %}
{%- endfor %}
{%- else %}
arguments: []
{%- endif %}
{%- if stdout %}
stdout: {{ stdout }}
{%- endif %}
description: |
{%- if help %}
  {{ help|indent(2) }}
{%- else %}
   TODO: Fill in description.
{%- endif %}
"""


def build(**kwds):
    """Build up a :func:`ToolDescription` from supplid arguments."""
    if kwds.get("cwl"):
        builder = _build_cwl
    else:
        builder = _build_galaxy
    return builder(**kwds)


def _build_cwl(**kwds):
    _handle_help(kwds)
    _handle_requirements(kwds)

    command_io = CommandIO(**kwds)
    render_kwds = {
        "cwl_version": "cwl:draft-3",
        "help": kwds.get("help", ""),
        "containers": kwds.get("containers", []),
        "id": kwds.get("id"),
        "label": kwds.get("name"),
    }
    render_kwds.update(command_io.cwl_properties())

    contents = _render(render_kwds, template_str=CWL_TEMPLATE)
    macro_contents = None
    test_files = []
    return ToolDescription(
        contents,
        macro_contents,
        test_files
    )


def _build_galaxy(**kwds):
    # Test case to build up from supplied inputs and outputs, ultimately
    # ignored unless kwds["test_case"] is truthy.

    _handle_help(kwds)

    # process raw cite urls
    cite_urls = kwds.get("cite_url", [])
    del kwds["cite_url"]
    citations = map(UrlCitation, cite_urls)
    kwds["bibtex_citations"] = citations

    # handle requirements and containers
    _handle_requirements(kwds)

    command_io = CommandIO(**kwds)
    kwds["inputs"] = command_io.inputs
    kwds["outputs"] = command_io.outputs
    kwds["command"] = command_io.cheetah_template

    test_case = command_io.test_case()

    # finally wrap up tests
    tests, test_files = _handle_tests(kwds, test_case)
    kwds["tests"] = tests

    # Render tool content from template.
    contents = _render(kwds)

    macro_contents = None
    if kwds["macros"]:
        macro_contents = _render(kwds, MACROS_TEMPLATE)

    return ToolDescription(contents, macro_contents, test_files)


class CommandIO(object):

    def __init__(self, **kwds):
        command = _find_command(kwds)
        cheetah_template = command

        # process raw inputs
        inputs = kwds.pop("input", [])
        inputs = list(map(Input, inputs or []))

        # alternatively process example inputs
        example_inputs = kwds.pop("example_input", [])
        for i, input_file in enumerate(example_inputs or []):
            name = "input%d" % (i + 1)
            inputs.append(Input(input_file, name=name, example=True))
            cheetah_template = _replace_file_in_command(cheetah_template, input_file, name)

        # handle raw outputs (from_work_dir ones) as well as named_outputs
        outputs = kwds.pop("output", [])
        outputs = list(map(Output, outputs or []))

        named_outputs = kwds.pop("named_output", [])
        for named_output in (named_outputs or []):
            outputs.append(Output(name=named_output, example=False))

        # handle example outputs
        example_outputs = kwds.pop("example_output", [])
        for i, output_file in enumerate(example_outputs or []):
            name = "output%d" % (i + 1)
            from_path = output_file
            use_from_path = True
            if output_file in cheetah_template:
                # Actually found the file in the command, assume it can
                # be specified directly and skip from_work_dir.
                use_from_path = False
            output = Output(name=name, from_path=from_path,
                            use_from_path=use_from_path, example=True)
            outputs.append(output)
            cheetah_template = _replace_file_in_command(cheetah_template, output_file, output.name)

        self.inputs = inputs
        self.outputs = outputs
        self.command = command
        self.cheetah_template = cheetah_template

    def example_input_names(self):
        for input in self.inputs:
            if input.example:
                yield input.input_description

    def example_output_names(self):
        for output in self.outputs:
            if output.example:
                yield output.example_path

    def cwl_lex_list(self):
        if not self.command:
            return []

        command_parts = shlex.split(self.command)
        parse_list = []

        input_count = 0
        output_count = 0

        index = 0

        prefixed_parts = []
        while index < len(command_parts):
            value = command_parts[index]
            eq_split = value.split("=")

            prefix = None
            if not _looks_like_start_of_prefix(index, command_parts):
                index += 1
            elif len(eq_split) == 2:
                prefix = Prefix(eq_split[0] + "=", False)
                value = eq_split[1]
                index += 1
            else:
                prefix = Prefix(value, True)
                value = command_parts[index + 1]
                index += 2
            prefixed_parts.append((prefix, value))

        for position, (prefix, value) in enumerate(prefixed_parts):
            if value in self.example_input_names():
                input_count += 1
                input = CwlInput("input%d" % input_count, position, prefix)
                parse_list.append(input)
            elif value in self.example_output_names():
                output_count += 1
                output = CwlOutput("output%d" % output_count, position, prefix)
                parse_list.append(output)
            elif prefix:
                param_id = prefix.prefix.lower().rstrip("=")
                type_ = param_type(value)
                input = CwlInput(param_id, position, prefix, type_=type_)
                parse_list.append(input)
            else:
                part = CwlCommandPart(value, position, prefix)
                parse_list.append(part)
        return parse_list

    def cwl_properties(self):
        base_command = []
        arguments = []
        inputs = []
        outputs = []

        lex_list = self.cwl_lex_list()

        index = 0
        while index < len(lex_list):
            token = lex_list[index]
            if isinstance(token, CwlCommandPart):
                base_command.append(token.value)
            else:
                break
            index += 1

        while index < len(lex_list):
            token = lex_list[index]
            if token.is_token(">"):
                break
            token.position = index - len(base_command) + 1
            if isinstance(token, CwlCommandPart):
                arguments.append(token)
            elif isinstance(token, CwlInput):
                inputs.append(token)
            elif isinstance(token, CwlOutput):
                token.glob = "$(inputs.%s)" % token.id
                outputs.append(token)

            index += 1

        stdout = None
        if index < len(lex_list):
            token = lex_list[index]
            if token.is_token(">") and (index + 1) < len(lex_list):
                output_token = lex_list[index + 1]
                if not isinstance(output_token, CwlOutput):
                    output_token = CwlOutput("std_out", None)

                output_token.glob = "out"
                output_token.require_filename = False
                outputs.append(output_token)
                stdout = "out"
                index += 2
            else:
                io.warn("Example command too complex, you will need to build it up manually.")

        return {
            "inputs": inputs,
            "outputs": outputs,
            "arguments": arguments,
            "base_command": base_command,
            "stdout": stdout,
        }

    def test_case(self):
        test_case = TestCase()
        for input in self.inputs:
            if input.example:
                test_case.params.append((input.name, input.input_description))

        for output in self.outputs:
            if output.example:
                test_case.outputs.append((output.name, output.from_path))

        return test_case


def _looks_like_start_of_prefix(index, parts):
    value = parts[index]
    if len(value.split("=")) == 2:
        return True
    if index + 1 == len(parts):
        return False
    next_value = parts[index + 1]
    next_value_is_not_start = (len(value.split("=")) != 2) and next_value[0] not in ["-", ">", "<", "|"]
    return value.startswith("-") and next_value_is_not_start


Prefix = namedtuple("Prefix", ["prefix", "separated"])


class CwlCommandPart(object):

    def __init__(self, value, position, prefix):
        self.value = value
        self.position = position
        self.prefix = prefix

    def is_token(self, value):
        return self.value == value


class CwlInput(object):

    def __init__(self, id, position, prefix, type_="File"):
        self.id = id
        self.position = position
        self.prefix = prefix
        self.type = type_

    def is_token(self, value):
        return False


class CwlOutput(object):

    def __init__(self, id, position, prefix):
        self.id = id
        self.position = position
        self.prefix = prefix
        self.glob = None
        self.require_filename = True

    def is_token(self, value):
        return False


def _render(kwds, template_str=TOOL_TEMPLATE):
    """ Apply supplied template variables to TOOL_TEMPLATE to generate
    the final tool.
    """
    return templates.render(template_str, **kwds)


def _replace_file_in_command(command, specified_file, name):
    """ Replace example file with cheetah variable name in supplied command
    or command template. Be sure to quote the name.
    """
    # TODO: check if the supplied variant was single quoted already.
    if '"%s"' % specified_file in command:
        # Sample command already wrapped filename in double quotes
        command = command.replace(specified_file, '$%s' % name)
    elif (" %s " % specified_file) in (" " + command + " "):
        # In case of spaces, best to wrap filename in double quotes
        command = command.replace(specified_file, '"$%s"' % name)
    else:
        command = command.replace(specified_file, '$%s' % name)
    return command


def _handle_help(kwds):
    """ Convert supplied help parameters into a help variable for template.
    If help_text is supplied, use as is. If help is specified from a command,
    run the command and use that help text.
    """
    help_text = kwds.get("help_text")
    if not help_text:
        help_from_command = kwds.get("help_from_command")
        if help_from_command:
            p = subprocess.Popen(
                help_from_command,
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT
            )
            help_text = p.communicate()[0]

    del kwds["help_text"]
    del kwds["help_from_command"]

    kwds["help"] = help_text


def _handle_tests(kwds, test_case):
    """ Given state built up from handling rest of arguments (test_case) and
    supplied kwds - build tests for template and corresponding test files.
    """
    test_files = []
    if kwds["test_case"]:
        tests = [test_case]
        test_files.extend(map(lambda x: x[1], test_case.params))
        test_files.extend(map(lambda x: x[1], test_case.outputs))
    else:
        tests = []
    return tests, test_files


def _handle_requirements(kwds):
    """ Convert requirements and containers specified from the command-line
    into abstract format for consumption by the template.
    """
    requirements = kwds["requirement"]
    del kwds["requirement"]
    requirements = map(Requirement, requirements or [])

    container = kwds["container"]
    del kwds["container"]
    containers = map(Container, container or [])

    kwds["requirements"] = requirements
    kwds["containers"] = containers


def _find_command(kwds):
    """Find base command from supplied arguments or just return None.

    If no such command was supplied (template will just replace this
    with a TODO item).
    """
    command = kwds.get("command")
    if not command:
        command = kwds.get("example_command", None)
        if command:
            del kwds["example_command"]
    return command


class UrlCitation(object):

    def __init__(self, url):
        self.url = url

    def __str__(self):
        if "github.com" in self.url:
            return self._github_str()
        else:
            return self._url_str()

    def _github_str(self):
        url = self.url
        title = url.split("/")[-1]
        return '''
@misc{github%s,
  author = {LastTODO, FirstTODO},
  year = {TODO},
  title = {%s},
  publisher = {GitHub},
  journal = {GitHub repository},
  url = {%s},
}''' % (title, title, url)

    def _url_str(self):
        url = self.url
        return '''
@misc{renameTODO,
  author = {LastTODO, FirstTODO},
  year = {TODO},
  title = {TODO},
  url = {%s},
}''' % (url)


class ToolDescription(object):
    """An description of the tool and related files to create."""

    def __init__(self, contents, macro_contents=None, test_files=[]):
        self.contents = contents
        self.macro_contents = macro_contents
        self.test_files = test_files


class Input(object):

    def __init__(self, input_description, name=None, example=False):
        parts = input_description.split(".")
        name = name or parts[0]
        if len(parts) > 0:
            datatype = ".".join(parts[1:])
        else:
            datatype = "data"

        self.input_description = input_description
        self.example = example
        self.name = name
        self.datatype = datatype

    def __str__(self):
        template = '<param type="data" name="{0}" format="{1}" />'
        return template.format(self.name, self.datatype)


class Output(object):

    def __init__(self, from_path=None, name=None, use_from_path=False, example=False):
        if from_path:
            parts = from_path.split(".")
            name = name or parts[0]
            if len(parts) > 1:
                datatype = ".".join(parts[1:])
            else:
                datatype = "data"
        else:
            name = name
            datatype = "data"

        self.name = name
        self.datatype = datatype
        if use_from_path:
            self.from_path = from_path
        else:
            self.from_path = None
        self.example = example
        if example:
            self.example_path = from_path

    def __str__(self):
        if self.from_path:
            return self._from_path_str()
        else:
            return self._named_str()

    def _from_path_str(self):
        template = '<data name="{0}" format="{1}" from_work_dir="{2}" />'
        return template.format(self.name, self.datatype, self.from_path)

    def _named_str(self):
        template = '<data name="{0}" format="{1}" />'
        return template.format(self.name, self.datatype)


class Requirement(object):

    def __init__(self, requirement):
        parts = requirement.split("@", 1)
        if len(parts) > 1:
            name = parts[0]
            version = "@".join(parts[1:])
        else:
            name = parts[0]
            version = None
        self.name = name
        self.version = version

    def __str__(self):
        base = '<requirement type="package"{0}>{1}</requirement>'
        if self.version is not None:
            attrs = ' version="{0}"'.format(self.version)
        else:
            attrs = ''
        return base.format(attrs, self.name)


def param_type(value):
    if re.match("^\d+$", value):
        return "int"
    elif re.match("^\d+?\.\d+?$", value):
        return "float"
    else:
        return "string"


class Container(object):

    def __init__(self, image_id):
        self.type = "docker"
        self.image_id = image_id

    def __str__(self):
        template = '<container type="{0}">{1}</container>'
        return template.format(self.type, self.image_id)


class TestCase(object):

    def __init__(self):
        self.params = []
        self.outputs = []


__all__ = [
    "build",
    "ToolDescription"
]
