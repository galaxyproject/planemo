import subprocess
from planemo import templates


TOOL_TEMPLATE = """<tool id="{{id}}" name="{{name}}" version="{{version}}">
{%- if description %}
    <description>{{ description }}</description>
{% endif %}
{%- if macros %}
    <macros>
        <import>macros.xml</import>
    </macros>
    <expand macro="requirements" />
    <expand macro="stdio" />
{%- if version_command %}
    <expand macro="version_command" />
{% endif %}
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
{% endif %}
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
{%- if doi %}
    <citations>
{%- for single_doi in doi %}
        <citation type="doi">{{ single_doi }}</citation>
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


def build(**kwds):
    # Test case to build up from supplied inputs and outputs, ultimately
    # ignored unless kwds["test_case"] is truthy.
    test_case = TestCase()

    command = _find_command(kwds)

    _handle_help(kwds)

    # process raw inputs
    inputs = kwds.get("input", [])
    del kwds["input"]
    inputs = list(map(Input, inputs or []))

    # alternatively process example inputs
    example_inputs = kwds["example_input"]
    del kwds["example_input"]
    for i, input_file in enumerate(example_inputs or []):
        name = "input%d" % (i + 1)
        inputs.append(Input(input_file, name=name))
        test_case.params.append((name, input_file))
        command = _replace_file_in_command(command, input_file, name)

    # handle raw outputs (from_work_dir ones) as well as named_outputs
    outputs = kwds.get("output", [])
    del kwds["output"]
    outputs = list(map(Output, outputs or []))

    named_outputs = kwds.get("named_output", [])
    del kwds["named_output"]
    for named_output in (named_outputs or []):
        outputs.append(Output(name=named_output))

    # handle example outputs
    example_outputs = kwds["example_output"]
    del kwds["example_output"]
    for i, output_file in enumerate(example_outputs or []):
        name = "output%d" % (i + 1)
        from_path = output_file
        if output_file in command:
            # Actually found the file in the command, assume it can
            # be specified directly and skip from_work_dir.
            from_path = None
        output = Output(name=name, from_path=from_path)
        outputs.append(output)
        test_case.outputs.append((name, output_file))
        command = _replace_file_in_command(command, output_file, output.name)

    kwds["inputs"] = inputs
    kwds["outputs"] = outputs

    # handle requirements and containers
    _handle_requirements(kwds)

    kwds["command"] = command

    # finally wrap up tests
    tests, test_files = _handle_tests(kwds, test_case)
    kwds["tests"] = tests

    # Render tool content from template.
    contents = _render(kwds)

    macro_contents = None
    if kwds["macros"]:
        macro_contents = _render(kwds, MACROS_TEMPLATE)

    return ToolDescription(contents, macro_contents, test_files)


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
    else:
        # In case of spaces, best to wrap filename in double quotes
        command = command.replace(specified_file, '"$%s"' % name)
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
    """ Find base command from supplied arguments or just return None if no
    such command was supplied (template will just replace this with TODO
    item).
    """
    command = kwds.get("command")
    if not command:
        command = kwds.get("example_command", None)
        if command:
            del kwds["example_command"]
    return command


class ToolDescription(object):

    def __init__(self, contents, macro_contents, test_files):
        self.contents = contents
        self.macro_contents = macro_contents
        self.test_files = test_files


class Input(object):

    def __init__(self, input_description, name=None):
        parts = input_description.split(".")
        name = name or parts[0]
        if len(parts) > 0:
            datatype = ".".join(parts[1:])
        else:
            datatype = "data"

        self.name = name
        self.datatype = datatype

    def __str__(self):
        template = '<param type="data" name="{0}" format="{1}" />'
        return template.format(self.name, self.datatype)


class Output(object):

    def __init__(self, from_path=None, name=None):
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
        self.from_path = from_path

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
