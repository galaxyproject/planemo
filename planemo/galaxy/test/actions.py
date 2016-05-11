"""Actions related to running and reporting on Galaxy-specific testing."""
import os

import click

from . import structures as test_structures
from planemo.io import error, info, warn, shell_join
from planemo.galaxy.run import (
    run_galaxy_command,
    setup_venv,
)
from planemo.reports import build_report


from galaxy.tools.deps.commands import shell

NO_XUNIT_REPORT_MESSAGE = ("Cannot locate xUnit report [%s] for tests - "
                           "required to build planemo report and summarize "
                           "tests.")
NO_JSON_REPORT_MESSAGE = ("Cannot locate JSON report [%s] for tests - "
                          "required to build planemo report and summarize "
                          "tests.")
REPORT_NOT_CHANGED = ("Galaxy failed to update test report [%s] for tests - "
                      "required to build planemo report and summarize "
                      "tests.")
NO_TESTS_MESSAGE = "No tests were executed - see Galaxy output for details."
ALL_TESTS_PASSED_MESSAGE = "All %d test(s) executed passed."
PROBLEM_COUNT_MESSAGE = ("There were problems with %d test(s) - out of %d "
                         "test(s) executed. See %s for detailed breakdown.")
GENERIC_PROBLEMS_MESSAGE = ("One or more tests failed. See %s for detailed "
                            "breakdown.")
GENERIC_TESTS_PASSED_MESSAGE = "No failing tests encountered."


def run_in_config(ctx, config, run=run_galaxy_command, **kwds):
    """Run Galaxy tests with the run_tests.sh command.

    The specified `config` object describes the context for tool
    execution.
    """
    config_directory = config.config_directory
    html_report_file = kwds["test_output"]

    job_output_files = kwds.get("job_output_files", None)
    if job_output_files is None:
        job_output_files = os.path.join(config_directory, "jobfiles")

    xunit_report_file = _xunit_state(kwds, config)
    xunit_report_file_tracker = _FileChangeTracker(xunit_report_file)
    structured_report_file = _structured_report_file(kwds, config)
    structured_report_file_tracker = _FileChangeTracker(structured_report_file)

    info("Testing using galaxy_root %s", config.galaxy_root)
    # TODO: Allow running dockerized Galaxy here instead.
    server_ini = os.path.join(config_directory, "galaxy.ini")
    config.env["GALAXY_CONFIG_FILE"] = server_ini
    config.env["GALAXY_TEST_VERBOSE_ERRORS"] = "true"
    config.env["GALAXY_TEST_SAVE"] = job_output_files

    cd_to_galaxy_command = "cd %s" % config.galaxy_root
    test_cmd = test_structures.GalaxyTestCommand(
        html_report_file,
        xunit_report_file,
        structured_report_file,
        failed=kwds.get("failed", False),
        installed=kwds.get("installed", False),
    ).build()
    setup_common_startup_args = ""
    if kwds.get("skip_venv", False):
        setup_common_startup_args = (
            'COMMON_STARTUP_ARGS=--skip-venv; '
            'export COMMON_STARTUP_ARGS; '
            'echo "Set COMMON_STARTUP_ARGS to ${COMMON_STARTUP_ARGS}"'
        )
    setup_venv_command = setup_venv(ctx, kwds)
    cmd = shell_join(
        cd_to_galaxy_command,
        setup_common_startup_args,
        setup_venv_command,
        test_cmd,
    )
    action = "Testing tools"
    return_code = run(
        ctx,
        cmd,
        config.env,
        action
    )
    if kwds.get('update_test_data', False):
        update_cp_args = (job_output_files, config.test_data_dir)
        shell('cp -r "%s"/* "%s"' % update_cp_args)

    _check_test_outputs(xunit_report_file_tracker, structured_report_file_tracker)

    test_results = test_structures.GalaxyTestResults(
        structured_report_file,
        xunit_report_file,
        html_report_file,
        return_code,
    )

    test_data = test_results.structured_data
    handle_reports(ctx, test_data, kwds)
    _handle_summary(
        test_results,
        **kwds
    )

    return return_code


def handle_reports(ctx, test_data, kwds):
    """Write reports based on user specified kwds."""
    exceptions = []
    for report_type in ["html", "markdown", "text"]:
        try:
            _handle_test_output_file(
                ctx, report_type, test_data, kwds
            )
        except Exception as e:
            exceptions.append(e)
            continue

    if len(exceptions) > 0:
        raise exceptions[0]


def _handle_test_output_file(ctx, report_type, test_data, kwds):
    kwd_name = "test_output"
    if report_type != "html":
        kwd_name = "test_output_%s" % report_type

    path = kwds.get(kwd_name, None)
    if path is None:
        message = "No file specified for %s, skipping test output." % kwd_name
        ctx.vlog(message)
        return

    try:
        contents = build_report.build_report(
            test_data, report_type=report_type
        )
    except Exception:
        message = "Problem producing report file %s for %s" % (
            path, kwd_name
        )
        ctx.vlog(message, exception=True)
        raise

    try:
        with open(path, 'w') as handle:
            handle.write(contents)
    except Exception:
        message = "Problem writing output file %s for %s" % (
            kwd_name, path
        )
        ctx.vlog(message, exception=True)
        raise


def _handle_summary(
    test_results,
    **kwds
):
    summary_style = kwds.get("summary")
    if summary_style == "none":
        return

    if test_results.has_details:
        _summarize_tests_full(
            test_results,
            **kwds
        )
    else:
        if test_results.exit_code:
            warn(GENERIC_PROBLEMS_MESSAGE % test_results.output_html_path)
        else:
            info(GENERIC_TESTS_PASSED_MESSAGE)


def _summarize_tests_full(
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
        _summarize_test_case(structured_data_tests, testcase_el, **kwds)


def passed(xunit_testcase_el):
    did_pass = True
    for child_el in list(xunit_testcase_el):
        if child_el.tag in ["failure", "error"]:
            did_pass = False
    return did_pass


def _summarize_test_case(structured_data, testcase_el, **kwds):
    summary_style = kwds.get("summary")
    test_id = test_structures.case_id(testcase_el)
    if not passed(testcase_el):
        state = click.style("failed", bold=True, fg='red')
    else:
        state = click.style("passed", bold=True, fg='green')
    click.echo(test_id.label + ": " + state)
    if summary_style != "minimal":
        _print_command_line(structured_data, test_id)


def _print_command_line(structured_data, test_id):
    try:
        test = [d for d in structured_data if d["id"] == test_id.id][0]["data"]
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


def _check_test_outputs(
    xunit_report_file_tracker,
    structured_report_file_tracker
):
    if not os.path.exists(xunit_report_file_tracker.path):
        message = NO_XUNIT_REPORT_MESSAGE % xunit_report_file_tracker.path
        error(message)
        raise Exception(message)

    if not os.path.exists(structured_report_file_tracker.path):
        message = NO_JSON_REPORT_MESSAGE % structured_report_file_tracker.path
        error(message)
        raise Exception(message)

    if not xunit_report_file_tracker.changed():
        message = REPORT_NOT_CHANGED % xunit_report_file_tracker.path
        error(message)
        raise Exception(message)

    if not structured_report_file_tracker.changed():
        message = REPORT_NOT_CHANGED % structured_report_file_tracker.path
        error(message)
        raise Exception(message)


def _xunit_state(kwds, config):
    # This has been supported in Galaxy for well over a year, just going to assume
    # it from here on out.
    # xunit_supported = True
    # if shell("grep -q xunit '%s'/run_tests.sh" % config.galaxy_root):
    #    xunit_supported = False

    xunit_report_file = kwds.get("test_output_xunit", None)
    if xunit_report_file is None:
        xunit_report_file = os.path.join(config.config_directory, "xunit.xml")

    return xunit_report_file


def _structured_report_file(kwds, config):
    # This has been supported in Galaxy for well over a year, just going to assume
    # it from here on out.
    # structured_data_supported = True
    # if shell("grep -q structured_data '%s'/run_tests.sh" % config.galaxy_root):
    #    structured_data_supported = False

    structured_report_file = None
    structured_report_file = kwds.get("test_output_json", None)
    if structured_report_file is None:
        conf_dir = config.config_directory
        structured_report_file = os.path.join(conf_dir, "structured_data.json")

    return structured_report_file


class _FileChangeTracker(object):

    def __init__(self, path):
        modification_time = None
        if os.path.exists(path):
            modification_time = os.path.getmtime(path)

        self.path = path
        self.modification_time = modification_time

    def changed(self):
        if self.modification_time:
            new_modification_time = os.path.getmtime(self.path)
            return self.modification_time != new_modification_time
        else:
            return os.path.exists(self.path)


__all__ = [
    "run_in_config",
    "handle_reports",
]
