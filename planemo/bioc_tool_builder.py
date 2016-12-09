"""Tool builder for R-bioc tools."""

import os

import yaml

from planemo.conda_recipes import write_bioconda_recipe

from planemo.tool_builder import (
    _find_command,
    _handle_help,
    _handle_tests,
    _render,
    _replace_file_in_command,
    Input,
    MACROS_TEMPLATE,
    TestCase,
    ToolDescription,
    UrlCitation,
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
    # inputs are in a list
    inputs = kwds.get("input", [])
    del kwds["input"]

    # TODO: DEPRICATE HANDLING OF EXAMPLE INPUT
    example_inputs = kwds["example_input"]
    del kwds["example_input"]

    rscript_data = kwds["rscript_data"]

    # Rscript inputs
    if bool(rscript_data):
        input_dict = rscript_data.get('inputs')  # dictionary of input parameters
        inputs = input_dict.values()[0]

    def edit_params(iodata):
        return iodata.split("/")[-1]
    inputs = map(edit_params, inputs)
    param_set = inputs
    inputs = list(map(Input, inputs or []))

    # TODO: DEPRICATE HANDLING OF EXAMPLE INPUT
    if not bool(rscript_data) and example_inputs:
        '''If no Rscript data is found but example_inputs are given - this should not happen'''
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
        outputs = output_dict.values()[0]
    # Add all parameters to a list called param_set
    outputs = map(edit_params, outputs)
    param_set.extend(outputs)

    outputs = list(map(Output, outputs or []))

    named_outputs = kwds.get("named_output", [])
    del kwds["named_output"]
    for named_output in (named_outputs or []):
        outputs.append(Output(name=named_output))

    # DEPRICATED HANDLING OF EXAMPLE OUTPUT
    # TODO: handle example outputs
    # if kwds.get("example_ouput"):
        # example_outputs = kwds["example_output"]
        # del kwds["example_output"]

    kwds['inputs'] = inputs
    kwds["outputs"] = outputs

    # handle requirements and containers
    if bool(rscript_data):
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

    # Edit command before sending it into the kwds dictionary
    command = _parse_command_rbioc(command, param_set)
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


def _parse_command_rbioc(command, param_set):
    """
    Parse command in Rbioc.

    Find a r-bioc command and replace the inputs and
    outputs with appropriate galaxy template.
    """
    cmd = command.split(" ")
    count = 0
    for i in range(len(cmd)):
        if "--" in cmd[i]:
            cmd[i + 1] = "$" + param_set[count].split(".")[0]
            count = count + 1
    return " ".join(cmd)


def _handle_requirements(kwds):
    """Handle Requirements in R-Bioc tools.

    Convert requirements specified from the command-line
    into abstract format for consumption by the template.
    """
    requirements = kwds["requirements"]
    # Handle Bioc/R requirements
    if kwds.get('rversion') and kwds.get("requirements"):
        rversion = kwds['rversion']
        rversion = rversion.replace(" ", "@")
        requirements.append(rversion)
        requirements.reverse()  # FIXME: Flips list so that R is the first in the requirements list
    bioconda_path = kwds["bioconda_path"]
    del kwds["requirements"]
    requirements = requirements or []

    requirements = [Requirement(req, bioconda_path=bioconda_path) for req in requirements]
    kwds["requirements"] = requirements


class Output(object):
    """Output class for R bioc tools."""

    def __init__(self, from_path=None, name=None, use_from_path=False):
        """Initialize Output class."""
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
        """Method override str function in Output class."""
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
    """Requirements class for R-Bioc tools."""

    def __init__(self, requirement, bioconda_path=None, update=False):
        """Initialize Requirements class for Bioc tool builder."""
        parts = requirement.split("@", 1)
        # Get version from requirements, if version not given
        if len(parts) > 1:
            name = parts[0]
            version = "@".join(parts[1:])
        else:
            name = parts[0]
            version = None
        if name == "R":
            self.name = name
            self.version = version
        else:
            # Write biconda recipe with given requirement
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
        """Format requirements string with correct version and name."""
        base = '<requirement type="package"{0}>{1}</requirement>'
        if self.version is not None:
            attrs = ' version="{0}"'.format(self.version)
        else:
            attrs = ''
        return base.format(attrs, self.name)
