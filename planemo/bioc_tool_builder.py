
import yaml
import os
from planemo.conda import write_bioconda_recipe
from planemo.tool_builder import (
    MACROS_TEMPLATE,
    ToolDescription,
    UrlCitation,
    Input,
    TestCase,
    _handle_help,
    _render,
    _replace_file_in_command,
)


def build(**kwds):
    """Build up a :func:`ToolDescription` from supplid arguments."""
    test_case = TestCase()

    command = _find_command(kwds)

    # process raw cite urls
    cite_urls = kwds.get("cite_url", [])
    del kwds["cite_url"]
    citations = map(UrlCitation, cite_urls)
    kwds["bibtex_citations"] = citations

    # process raw inputs
    inputs = kwds.get("input", [])
    del kwds["input"]
    # alternatively process example inputs
    example_inputs = kwds["example_input"]
    del kwds["example_input"]

    # Rscript inputs
    rscript_data = kwds["rscript_data"]
    if bool(rscript_data):
        input_dict = rscript_data.get('inputs')
        inputs = list(input_dict.values())
    # print(inputs)
    inputs = list(map(Input, inputs or []))

    if not bool(rscript_data) and example_inputs:
        for i, input_file in enumerate(example_inputs or []):
            name = "input%d" % (i + 1)
            inputs.append(Input(input_file, name=name))
            test_case.params.append((name, input_file))
            command = _replace_file_in_command(command, input_file, name)

    # handle raw outputs (from_work_dir ones) as well as named_outputs
    outputs = kwds.get("output", [])
    del kwds["output"]

    if bool(rscript_data):
        output_dict = rscript_data.get('outputs')
        outputs = list(output_dict.values())
    # print(outputs)

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
        use_from_path = True
        if output_file in command:
            # Actually found the file in the command, assume it can
            # be specified directly and skip from_work_dir.
            use_from_path = False
        output = Output(name=name, from_path=from_path,
                        use_from_path=use_from_path)
        outputs.append(output)
        test_case.outputs.append((name, output_file))
        command = _replace_file_in_command(command, output_file, output.name)

    kwds["inputs"] = inputs
    kwds["outputs"] = outputs

    # handle requirements and containers
    if bool(rscript_data):
        # print(rscript_data.get('library'))
        if not type([rscript_data.get('library')]) is list:
            kwds['requirements'] = [rscript_data.get('library')]
        else:
            kwds['requirements'] = rscript_data.get('library')

    _handle_requirements(kwds)

    # Add help from requirements
    if not bool(rscript_data):
        req = kwds['requirements'][0]
        command_help = req.package_help + "\n \n" + req.package_url
        kwds['help_text'] = command_help

    # Handle help
    _handle_help(kwds)

    kwds["command"] = command

    # finally wrap up tests
    tests, test_files = _handle_tests(kwds, test_case)
    kwds["tests"] = tests

    # Render tool content from template.
    contents = _render(kwds)
    macro_contents = None
    if kwds["macros"]:
        macro_contents = _render(kwds, MACROS_TEMPLATE)

    return ToolDescription(
        contents,
        macro_contents,
        test_files=test_files
    )


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
    requirements = kwds["requirements"]
    bioconda_path = kwds["bioconda_path"]
    del kwds["requirements"]
    requirements = requirements or []

    requirements = [Requirement(req, bioconda_path=bioconda_path) for req in requirements]

    # container = kwds["container"]
    # del kwds["container"]
    # containers = map(Container, container or [])

    kwds["requirements"] = requirements
    # kwds["containers"] = containers


def _find_command(kwds):
    """ Find base command from supplied arguments or just return None if no
    such k was supplied (template will just replace this with TODO
    item).
    """
    command = kwds.get("command")
    if not command:
        command = kwds.get("example_command", None)
        if command:
            del kwds["example_command"]
    return command


class Output(object):

    def __init__(self, from_path=None, name=None, use_from_path=False):
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

    def __init__(self, requirement, bioconda_path=None, update=False):

        parts = requirement.split("@", 1)
        # Get version from requirements, if version not given
        if len(parts) > 1:
            name = parts[0]
            version = "@".join(parts[1:])
        else:
            name = parts[0]
            version = None
        # Write biconda recipe with give requirement
        if bioconda_path is None:
            bioconda_path = os.path.expanduser("~")
        write_bioconda_recipe(name, True, update, bioconda_path)

        recipe_path = os.path.join(bioconda_path,
                                   "bioconda-recipes",
                                   "recipes",
                                   "bioconductor-" + name.lower(),
                                   "meta.yaml")
        if not os.path.exists(recipe_path):
            recipe_path = os.path.join(bioconda_path,
                                       "bioconda-recipes",
                                       "recipes",
                                       "r-" + name.lower(),
                                       "meta.yaml")
        with open(recipe_path, 'r') as f:
            doc = yaml.load(f)
            if not version:
                version = doc["package"]["version"]
            package_url = doc["about"]["home"]
            package_help = doc["about"]["summary"]
        self.name = name
        self.version = version
        self.package_url = package_url
        self.package_help = package_help

    def __str__(self):
        base = '<requirement type="package"{0}>{1}</requirement>'
        if self.version is not None:
            attrs = ' version="{0}"'.format(self.version)
        else:
            attrs = ''
        return base.format(attrs, self.name)
