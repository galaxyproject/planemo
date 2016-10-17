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
                    "(e.g. --example_command 'Rscript /full/path/to/my_r_tool.R --input input.csv --output output.csv'). "
                    "This option cannot be used with --command. Instead, this option "
                    "should be used --example_input and --example_output.")
REQUIREMENT_HELP = ("Name of the bioconductor package. "
                    "Requirements will be set using bioconda. (e.g. --requirement 'motifbreakR') ")


@click.command("bioc_tool_init")
@options.tool_init_id_option(prompt=False)
@options.force_option(what="tool")
@options.tool_init_tool_option()
@options.tool_init_name_option(
    prompt=False,
    help=NAME_OPTION_HELP,
)
@options.tool_init_version_option()
@options.tool_init_description_option()
@options.tool_init_command_option()
@options.tool_init_example_command_option(help=EXAMPLE_CMD_HELP)
@options.tool_init_requirement_option(help=REQUIREMENT_HELP)
@options.tool_init_example_input_option()
@options.tool_init_example_output_option()
@options.tool_init_named_output_option()
@options.tool_init_input_option()
@options.tool_init_output_option()
@options.tool_init_help_text_option()
@options.tool_init_help_from_command_option()
@options.tool_init_doi_option()
@options.tool_init_cite_url_option()
@options.tool_init_test_case_option()
@options.tool_init_macros_option()
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
@click.option(
    "--rscript",
    type=click.Path(exists=True),
    default=None,
    prompt=False,
    help=("Name of an R script - designed per Galaxy R tool "
          "best practices - from which to create a tool definition file."
          " (e.g. --rscript 'file.R') ")
)
@command_function
def cli(ctx, **kwds):
    """Generate a bioconductor tool outline from supplied arguments."""
    invalid = _validate_kwds(kwds)
    if kwds.get("rscript") and kwds.get("example_command"):
        rscript = kwds["rscript"]
        example_command = kwds["example_command"]
        rscript_data = rscript_parse.parse_rscript(rscript, example_command)
        # print(rscript_data)
        # Get name replace .R, or .r
        kwds['name'] = kwds.get("name")
        kwds['id'] = rscript.split("/")[-1].replace(".R", "")
        kwds['rscript_data'] = rscript_data
    else:  # if no rscript
        info("No Rscript found. Must provide correct planemo arguments.")

    if invalid:
        ctx.exit(invalid)

    tool_description = bioc_tool_builder.build(**kwds)
    tool_builder.write_tool_description(
        ctx, tool_description, **kwds
    )


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

    if not_exclusive("help_text", "help_from_command"):
        return 1
    if not_exclusive("command", "example_command"):
        return 1
    if not_exclusive("rscript", "requirement"):
        return 1
    if not_specifing_dependent_option("example_input", "example_command"):
        return 1
    if not_specifing_dependent_option("example_output", "example_command"):
        return 1
    if not_specifing_dependent_option("test_case", "example_command"):
        return 1
    if not_specifing_dependent_option("rscript", "example_command"):
        return 1
    return 0
