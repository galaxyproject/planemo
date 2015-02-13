import os
import sys

import click

from planemo.cli import pass_context
from planemo.io import info, warn
from planemo import options
from planemo import galaxy_config
from planemo import galaxy_run
from planemo import galaxy_test
from planemo.reports import build_report


from galaxy.tools.deps.commands import shell

XUNIT_UPGRADE_MESSAGE = ("This version of Galaxy does not support xUnit - "
                         "please update to newest development brach.")
NO_XUNIT_MESSAGE = ("Cannot locate xUnit report option for tests - update "
                    "Galaxy for more detailed breakdown.")
NO_JSON_MESSAGE = ("Cannot locate json report option for tests - update "
                   "Galaxy for more detailed breakdown.")
NO_TESTS_MESSAGE = "No tests were executed - see Galaxy output for details."
ALL_TESTS_PASSED_MESSAGE = "All %d test(s) executed passed."
PROBLEM_COUNT_MESSAGE = ("There were problems with %d test(s) - out of %d "
                         "test(s) executed. See %s for detailed breakdown.")
GENERIC_PROBLEMS_MESSAGE = ("One or more tests failed. See %s for detailed "
                            "breakdown.")
GENERIC_TESTS_PASSED_MESSAGE = "No failing tests encountered."

RUN_TESTS_CMD = (
    "sh run_tests.sh --report_file %s %s %s "
    " functional.test_toolbox"
)

OUTPUT_DFEAULTS = {
    "output": "tool_test_output.html",
    "output_json": "tool_test_output.json",
    "output_xunit": None,
}


@click.command('test')
@options.optional_tools_arg()
@click.option(
    "--test_output",
    type=click.Path(file_okay=True, resolve_path=True),
    help=("Output test report (HTML - for humans) defaults to "
          "tool_test_output.html."),
    default=None,
)
@click.option(
    "--test_output_xunit",
    type=click.Path(file_okay=True, resolve_path=True),
    help="Output test report (xUnit style - for computers).",
    default=None,
)
@click.option(
    "--test_output_json",
    type=click.Path(file_okay=True, resolve_path=True),
    help=("Output test report (planemo json) defaults to "
          "tool_test_output.json."),
    default=None,
)
@click.option(
    "--job_output_files",
    type=click.Path(file_okay=False, resolve_path=True),
    help="Write job outputs to specified directory.",
    default=None,
)
@click.option(
    "--update_test_data",
    is_flag=True,
    help="Update test-data directory with job outputs (normally written to "
         "directory --job_output_files if specified.)"
)
@click.option(
    "--summary",
    type=click.Choice(['none', 'minimal', 'compact']),
    default="minimal",
    help=("Summary style printed to planemo's standard output (see output "
          "reports for more complete summary). Set to 'none' to disable "
          "completely.")
)
@options.galaxy_root_option()
@options.install_galaxy_option()
@options.test_data_option()
@options.tool_data_table_option()
@options.dependency_resolvers_option()
@options.job_config_option()
@options.tool_dependency_dir_option()
@options.brew_dependency_resolution()
@options.shed_dependency_resolution()
@pass_context
def cli(ctx, path, **kwds):
    """Run the tests in the specified tool tests in a Galaxy instance.

    All referenced tools (by default all the tools in the current working
    directory) will be tested and the results quickly summarized.

    To run these tests planemo needs a Galaxy instance to utilize, planemo
    will search parent directories to see if any is a Galaxy instance
    - but one can pick the Galaxy instance to use with the --galaxy_root
    option or force planemo to download a disposable instance with the
    ``--install_galaxy`` flag.

    In additon to to quick summary printed to the console - various detailed
    output summaries can be configured. ``tool_test_output.html`` (settable
    via ``--test_output``) will contain a human consumable HTML report
    describing the test run. A JSON file (settable via ``--test_output_json``
    and defaulting to ``tool_test_output.json``) will also be created. These
    files can can be disabled by passing in empty arguments or globally by
    setting the values ``default_test_output`` and/or
    ``default_test_output_json`` in ``~/.planemo.yml`` to ``null``. For
    continuous integration testing a xUnit-style report can be confiured using
    the ``--test_output_xunit``.

    planemo uses temporarily generated config files and environment variables
    to attempt to shield this execution of Galaxy from manually launched runs
    against that same Galaxy root - but this may not be bullet proof yet so
    please careful and do not try this against production Galaxy instances.
    """
    for name, default in OUTPUT_DFEAULTS.items():
        __populate_default_output(ctx, name, kwds, default)

    kwds["for_tests"] = True
    with galaxy_config.galaxy_config(ctx, path, **kwds) as config:
        config_directory = config.config_directory
        html_report_file = kwds["test_output"]

        job_output_files = kwds.get("job_output_files", None)
        if job_output_files is None:
            job_output_files = os.path.join(config_directory, "jobfiles")

        xunit_supported, xunit_report_file = __xunit_state(kwds, config)
        structured_report_file = __structured_report_file(kwds, config)

        info("Testing using galaxy_root %s", config.galaxy_root)
        # TODO: Allow running dockerized Galaxy here instead.
        server_ini = os.path.join(config_directory, "galaxy.ini")
        config.env["GALAXY_CONFIG_FILE"] = server_ini
        config.env["GALAXY_TEST_VERBOSE_ERRORS"] = "true"
        config.env["GALAXY_TEST_SAVE"] = job_output_files

        cd_to_galaxy_command = "cd %s" % config.galaxy_root
        cmd = "; ".join([
            cd_to_galaxy_command,
            galaxy_run.ACTIVATE_COMMAND,  # TODO: this should be moved to
                                          # run_tests.sh to match run.sh.
            __run_tests_cmd(
                html_report_file,
                xunit_report_file,
                structured_report_file,
            ),
        ])
        action = "Testing tools"
        return_code = galaxy_run.run_galaxy_command(
            ctx,
            cmd,
            config.env,
            action
        )
        if kwds.get('update_test_data', False):
            update_cp_args = (job_output_files, config.test_data_dir)
            shell('cp -r "%s"/* "%s"' % update_cp_args)

        if xunit_report_file and (not os.path.exists(xunit_report_file)):
            warn(NO_XUNIT_MESSAGE)
            xunit_report_file = None

        test_results = galaxy_test.GalaxyTestResults(
            structured_report_file,
            xunit_report_file,
            html_report_file,
            return_code,
        )

        try:
            test_data = test_results.structured_data
            new_report = build_report.build_report(test_data)
            open(test_results.output_html_path, "w").write(new_report)
        except Exception:
            pass

        __handle_summary(
            test_results,
            **kwds
        )

        if return_code:
            sys.exit(1)


def __populate_default_output(ctx, type, kwds, default):
    kwd_key = "test_%s" % type
    kwd_value = kwds.get(kwd_key, None)
    if kwd_value is None:
        global_config = ctx.global_config
        global_config_key = "default_test_%s" % type
        if global_config_key in global_config:
            default_value = global_config[global_config_key]
        else:
            default_value = default

        if default_value:
            default_value = os.path.abspath(default_value)
        kwds[kwd_key] = default_value


def __handle_summary(
    test_results,
    **kwds
):
    summary_style = kwds.get("summary")
    if summary_style == "none":
        return

    if test_results.has_details:
        __summarize_tests_full(
            test_results,
            **kwds
        )
    else:
        if test_results.exit_code:
            warn(GENERIC_PROBLEMS_MESSAGE % test_results.output_html_path)
        else:
            info(GENERIC_TESTS_PASSED_MESSAGE)


def __summarize_tests_full(
    test_results,
    **kwds
):
    num_tests = test_results.num_tests
    num_problems = test_results.num_problems

    if num_tests == 0:
        warn(NO_TESTS_MESSAGE)
        return

    if num_problems == 0:
        info(ALL_TESTS_PASSED_MESSAGE % num_tests)

    if num_problems:
        html_report_file = test_results.output_html_path
        message_args = (num_problems, num_tests, html_report_file)
        message = PROBLEM_COUNT_MESSAGE % message_args
        warn(message)

    for testcase_el in test_results.xunit_testcase_elements:
        structured_data_tests = test_results.structured_data_tests
        __summarize_test_case(structured_data_tests, testcase_el, **kwds)


def __summarize_test_case(structured_data, testcase_el, **kwds):
    summary_style = kwds.get("summary")
    name_raw = testcase_el.attrib["name"]
    tool_and_num = name_raw.split("TestForTool_", 1)[-1]
    if "test_tool_" in tool_and_num:
        tool, num = tool_and_num.split(".test_tool_", 1)
        try:
            num = int(num)
        except ValueError:
            pass
        # Tempted to but something human friendly in here like
        # num + 1 - but then it doesn't match HTML report.
        label = "%s[%s]" % (tool, num)
    else:
        label = tool_and_num
    passed = len(list(testcase_el)) == 0
    if not passed:
        state = click.style("failed", bold=True, fg='red')
    else:
        state = click.style("passed", bold=True, fg='green')
    click.echo(label + ": " + state)
    if summary_style != "minimal":
        __print_command_line(structured_data, name_raw)


def __print_command_line(structured_data, test_id):
    try:
        test = [d for d in structured_data if d["id"] == test_id][0]["data"]
    except (KeyError, IndexError):
        # Failed to find structured data for this test - likely targetting
        # and older Galaxy version.
        return

    execution_problem = test.get("execution_problem", None)
    if execution_problem:
        click.echo("| command: *could not execute job, no command generated* ")
        return

    try:
        command = test["job"]["command_line"]
    except (KeyError, IndexError):
        click.echo("| command: *failed to determine command for job* ")
        return

    click.echo("| command: %s" % command)


def __run_tests_cmd(html_report_file, xunit_report_file, sd_report_file):
    if xunit_report_file:
        xunit_arg = "--xunit_report_file %s" % xunit_report_file
    else:
        xunit_arg = ""
    if sd_report_file:
        sd_arg = "--structured_data_report_file %s" % sd_report_file
    else:
        sd_arg = ""
    return RUN_TESTS_CMD % (html_report_file, xunit_arg, sd_arg)


def __xunit_state(kwds, config):
    xunit_supported = True
    if shell("grep -q xunit '%s'/run_tests.sh" % config.galaxy_root):
        xunit_supported = False

    xunit_report_file = kwds.get("test_output_xunit", None)
    if xunit_report_file is None and xunit_supported:
        xunit_report_file = os.path.join(config.config_directory, "xunit.xml")
    elif xunit_report_file is not None and not xunit_supported:
        warn(XUNIT_UPGRADE_MESSAGE)
        xunit_report_file = None

    return xunit_supported, xunit_report_file


def __structured_report_file(kwds, config):
    structured_data_supported = True
    if shell("grep -q structured_data '%s'/run_tests.sh" % config.galaxy_root):
        structured_data_supported = False

    structured_report_file = None
    structured_report_file = kwds.get("test_output_json", None)
    if structured_report_file is None and structured_data_supported:
        conf_dir = config.config_directory
        structured_report_file = os.path.join(conf_dir, "structured_data.json")
    elif structured_report_file is not None and not structured_data_supported:
        warn(NO_JSON_MESSAGE)
        structured_report_file = None

    return structured_report_file
