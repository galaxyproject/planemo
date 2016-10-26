"""Module describing the planemo ``bioc_tool_init`` command."""

import click

from planemo import bioc_tool_builder
# from planemo import bioconductor_skeleton
from planemo import io
from planemo import options
from planemo import rscript_parse
from planemo import tool_builder
from planemo.cli import command_function
from planemo.io import info

NAME_OPTION_HELP = "Name for new R/Bioconductor tool (user facing)."
EXAMPLE_CMD_HELP = ("Example command with paths to build Cheetah template from. "
                    "(e.g. --example_command 'Rscript /path/to/my_r_tool.R --input input.csv --output output.csv'). "
                    "This option cannot be used with --command. Instead, this option "
                    "should be used with --example_input and --example_output.")
REQUIREMENT_HELP = ("Name of the R/Bioconductor package. "
                    "Requirements will be set using Bioconda. (e.g. --requirement 'affy') ")


@click.command("bioc_tool_init")
@options.tool_init_id_option(prompt=False)
@options.force_option(what="tool")
@options.tool_init_tool_option()
@options.tool_init_name_option(
    prompt=False,
    help=NAME_OPTION_HELP,
)
@options.tool_init_description_option()
@options.tool_init_command_option()
@click.option(
    "--rscript",
    type=click.Path(exists=True),
    default=None,
    prompt=False,
    help=("Name of an R script from which to create a Tool definition file. "
          "Requires use of --input and --output arguments. "
          "(e.g. --rscript 'file.R') ")
)
@click.option(
    "--rversion",
    default="R 3.2.1",
    prompt=False,
    help=("R version this tool requries, if not given, "
          "the tool defaults to R 3.2.1, "
          "(eg: --rversion 'R 3.2.1'). This option adds the R dependency "
          "in the tool requirements.")
)
@options.tool_init_input_option()
@options.tool_init_output_option()
@options.tool_init_requirement_option(help=REQUIREMENT_HELP)
@options.tool_init_help_text_option()
@options.tool_init_doi_option()
@options.tool_init_cite_url_option()
@options.tool_init_version_option()
# TODO: The following are not actually used, remove later
@options.tool_init_help_from_command_option()
@options.tool_init_test_case_option()
@options.tool_init_macros_option()
@options.tool_init_named_output_option()
# Shares all basic Galaxy tool options with tool_init except version_command
# and container.
@click.option(
    "--bioconda_path",
    type=click.STRING,
    default=None,
    prompt=False,
    help=("Path to bioconda repository. "
          "If left empty, path will be made in home directory.")
)
# TODO: THESE FEATURES ARE NOT YET IMPLEMENTED
# @options.tool_init_example_command_option()
# @options.tool_init_example_input_option()
# @options.tool_init_example_output_option()
@command_function
def cli(ctx, **kwds):
    """Generate a bioconductor tool outline from supplied arguments."""
    invalid = _validate_kwds(kwds)
    _check_r_version(kwds)
    if kwds.get("command"):
        command = kwds["command"]
        rscript = command.split()[1]  # Name of Custom R file

    elif kwds.get("rscript") and kwds.get("input") and kwds.get("output"):
        rscript = kwds["rscript"]
        command = 'Rscript %s ' % rscript  # Build command from --rscript, --input, --output
        for i in kwds["input"]:
            command += '--input %s ' % i
        for o in kwds["output"]:
            command += '--output %s ' % o
    else:  # No --rscript/input/output and no --command given
        info("Need to supply EITHER a full command (--command) OR an R script (--rscript), input(s) (--input), and output(s) (--output).")
        ctx.exit(1)

    if invalid:
        ctx.exit(invalid)

    rscript_data = rscript_parse.parse_rscript(rscript, command)
    kwds['rscript_data'] = rscript_data
    kwds['rscript'] = rscript
    kwds['command'] = command
    kwds['name'] = kwds.get("name")
    kwds['id'] = rscript.split("/")[-1].replace(".R", "")  # Default: name of R script w/o extension

    # Assign input/output to kwds if --input/--output not used
    if not kwds['input']:
        new_inputs = ()
        for i in kwds['rscript_data']['inputs']['input']:
            new_inputs += (i,)
        kwds['input'] = new_inputs

    if not kwds['output']:
        new_outputs = ()
        for i in kwds['rscript_data']['outputs']['output']:
            new_outputs += (i,)
        kwds['output'] = new_outputs

    input_dict = rscript_data.get('inputs')
    inputs = list(input_dict.values())
    kwds['inputs'] = inputs

    # Set example_input/output to input/output for now
    # Probably can remove example_input/output in future
    kwds['example_input'] = kwds['input']
    kwds['example_output'] = kwds['output']

    # Build Tool definition file
    tool_description = bioc_tool_builder.build(**kwds)
    tool_builder.write_tool_description(
        ctx, tool_description, **kwds
    )


def _check_r_version(kwds):
    """Check if R version is supplied."""
    if kwds.get("rversion"):
        rversion = kwds["rversion"]
        print("R version {0} is being used.".format(rversion))
    elif kwds.get("rscript") and not kwds.get("rversion"):
        print("Default R 3.2.1 version will be put in the requirements")
    else:
        print("No R requirement given")


def _validate_kwds(kwds):
    def not_exclusive(x, y):
        if kwds.get(x) and kwds.get(y):
            io.error("Can only specify one of --%s and --%s" % (x, y))
            return True

    def not_specifing_dependent_option(x, y):
        if kwds.get(x) and not kwds.get(y):
            template = "Can only use the --%s option if also specifying --%s"
            message = template % (x, y)
            io.error(message)
            return True

    if not_specifing_dependent_option("rscript", "input"):
        return 1
    if not_specifing_dependent_option("rscript", "output"):
        return 1
    if not_exclusive("help_text", "help_from_command"):
        return 1
    if not_exclusive("command", "example_command"):
        return 1
    if not_exclusive("rscript", "requirement"):
        return 1
    if not_specifing_dependent_option("test_case", "example_command"):
        return 1
    return 0
