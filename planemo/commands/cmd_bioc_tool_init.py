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


# --input_format
# --output_format
# --advanced_options
@click.command("bioc_tool_init")
@options.tool_init_id_option(prompt=False)
@options.force_option(what="tool")
@options.tool_init_tool_option()
@options.tool_init_name_option(
    prompt=False,
    help="Name for new Bioconductor tool (user facing)",
)
@options.tool_init_version_option()
@options.tool_init_description_option()
@options.tool_init_command_option()
# TODO: Change this
@click.option(
    "--example_command",
    type=click.STRING,
    default=None,
    prompt=False,
    help=("Example to command with paths to build Cheetah template from "
          "(e.g. 'Rscript my_r_tool.R --input input.csv --output output.csv')"
          ". Option cannot be used with --command,"
          "should be used --example_input and --example_output."),
)
@click.option(
    "--requirement",
    type=click.STRING,
    default=None,
    multiple=True,
    prompt=False,
    help=("Give the name of the bioconductor package,"
          "requirements will be set using bioconda. eg: 'motifbreakR' ")
)
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
@click.option(
    "--bioconda_path",
    type=click.STRING,
    default=None,
    prompt=False,
    help=("Give path to bioconda repository,"
          " if left empty, path will be made in home directory")
)
@click.option(
    "--rscript",
    type=click.Path(exists=True),
    default=None,
    prompt=False,
    help=("Give an R Script, designed as per Galaxy R tool"
          "best practices, and create a tool definition file."
          "eg: planemo bioc_tool_init --rscript 'file.R' ")
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
        kwds['name'] = rscript.split("/")[-1].replace(".R", "")
        kwds['id'] = kwds.get("name")
        kwds['rscript_data'] = rscript_data
    else:  # if no rscript
        info("No Rscript found, must provide correct planemo arguments.")

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
