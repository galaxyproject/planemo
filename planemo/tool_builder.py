import subprocess
try:
    from jinja2 import Template
except ImportError:
    Template = None

NO_JINJA2_MESSAGE = ("This functionality requires Jinja2 but this library is "
                     "unavailable. Install with `pip install jinja2`.")
TOOL_TEMPLATE = """<tool id="{{id}}" name="{{name}}" version="{{version}}">
{%- if description %}
    <description>{{ description }}</description>
{% endif %}
    <stdio>
        <exit_code range="1:" />
    </stdio>
    <requirements>
{%- for requirement in requirements %}
        {{ requirement }}
{%- endfor %}
{%- for container in containers %}
        {{ container }}
{%- endfor %}
    </requirements>
{%- if version_command %}
    <version_command>{{ version_command }}</version_command>
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
</tool>
"""


def build(**kwds):
    if Template is None:
        raise Exception(NO_JINJA2_MESSAGE)

    command = kwds.get("command")
    if not command:
        command = kwds.get("example_command", None)
        if command:
            del kwds["example_command"]

    help_text = kwds.get("help_text")
    test_case = TestCase()

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

    # process raw inputs
    inputs = kwds.get("input", [])
    del kwds["input"]
    inputs = map(Input, inputs or [])

    example_inputs = kwds["example_input"]
    del kwds["example_input"]
    for i, input_file in enumerate(example_inputs or []):
        name = "input%d" % (i + 1)
        inputs.append(Input(input_file, name=name))
        test_case.params.append((name, input_file))
        command = command.replace(input_file, "$%s" % name)

    outputs = kwds.get("output", [])
    del kwds["output"]
    outputs = map(Output, outputs or [])

    named_outputs = kwds.get("named_output", [])
    del kwds["named_output"]
    for named_output in enumerate(named_outputs or []):
        outputs.append(Output(name=named_output))

    example_outputs = kwds["example_output"]
    del kwds["example_output"]
    for i, output_file in enumerate(example_outputs or []):
        name = "output%d" % (i + 1)
        output = Output(name=name, from_path=output_file)
        outputs.append(output)
        test_case.outputs.append((name, output_file))
        command = command.replace(output_file, "$%s" % output.name)

    requirements = kwds["requirement"]
    del kwds["requirement"]
    requirements = map(Requirement, requirements or [])

    container = kwds["container"]
    del kwds["container"]
    containers = map(Container, container or [])

    kwds["inputs"] = inputs
    kwds["outputs"] = outputs
    kwds["requirements"] = requirements
    kwds["containers"] = containers
    kwds["command"] = command
    test_files = []
    if kwds["test_case"]:
        kwds["tests"] = [test_case]
        test_files.extend(map(lambda x: x[1], test_case.params))
        test_files.extend(map(lambda x: x[1], test_case.outputs))
    else:
        kwds["tests"] = []
    template = Template(TOOL_TEMPLATE)
    contents = template.render(**kwds)

    return ToolDescription(contents, test_files)


class ToolDescription(object):

    def __init__(self, contents, test_files):
        self.contents = contents
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
