"""Actions related to running and reporting on Galaxy-specific testing."""
import json
import os

import click

from . import structures as test_structures
from planemo.test.results import get_dict_value
from planemo.io import error, info, warn, shell_join
from planemo.exit_codes import (
    EXIT_CODE_OK,
    EXIT_CODE_GENERIC_FAILURE,
    EXIT_CODE_NO_SUCH_TARGET,
)
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

    structured_data = test_results.structured_data
    return handle_reports_and_summary(
        ctx,
        structured_data,
        exit_code=test_results.exit_code,
        kwds=kwds
    )


def handle_reports_and_summary(ctx, structured_data, exit_code=None, kwds={}):
    """Produce reports and print summary, return 0 if tests passed.

    If ``exit_code`` is set - use underlying test source for return
    code and test success determination, otherwise infer from supplied
    test data.
    """
    handle_reports(ctx, structured_data, kwds)
    summary_exit_code = _handle_summary(
        structured_data,
        **kwds
    )
    return exit_code if exit_code is not None else summary_exit_code


def handle_reports(ctx, structured_data, kwds):
    """Write reports based on user specified kwds."""
    exceptions = []
    structured_report_file = kwds.get("test_output_json", None)
    if structured_report_file and not os.path.exists(structured_report_file):
        try:
            with open(structured_report_file, "w") as f:
                json.dump(structured_report_file, f)
        except Exception as e:
            exceptions.append(e)

    for report_type in ["html", "markdown", "text"]:
        try:
            _handle_test_output_file(
                ctx, report_type, structured_data, kwds
            )
        except Exception as e:
            exceptions.append(e)

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
    structured_data,
    **kwds
):
    summary_dict = get_dict_value("summary", structured_data)
    num_tests = get_dict_value("num_tests", summary_dict)
    num_failures = get_dict_value("num_failures", summary_dict)
    num_errors = get_dict_value("num_errors", summary_dict)
    num_problems = num_failures + num_errors

    summary_exit_code = EXIT_CODE_OK
    if num_problems > 0:
        summary_exit_code = EXIT_CODE_GENERIC_FAILURE
    elif num_tests == 0:
        summary_exit_code = EXIT_CODE_NO_SUCH_TARGET

    summary_style = kwds.get("summary")
    if summary_style != "none":
        if num_tests == 0:
            warn(NO_TESTS_MESSAGE)
        elif num_problems == 0:
            info(ALL_TESTS_PASSED_MESSAGE % num_tests)
        elif num_problems:
            html_report_file = kwds.get("test_output")
            message_args = (num_problems, num_tests, html_report_file)
            message = PROBLEM_COUNT_MESSAGE % message_args
            warn(message)

        _summarize_tests_full(
            structured_data,
            **kwds
        )

    return summary_exit_code


def _summarize_tests_full(
    structured_data,
    **kwds
):
    tests = get_dict_value("tests", structured_data)
    for test_case_data in tests:
        _summarize_test_case(test_case_data, **kwds)


def passed(xunit_testcase_el):
    did_pass = True
    for child_el in list(xunit_testcase_el):
        if child_el.tag in ["failure", "error"]:
            did_pass = False
    return did_pass


def _summarize_test_case(structured_data, **kwds):
    summary_style = kwds.get("summary")
    test_id = test_structures.case_id(
        raw_id=get_dict_value("id", structured_data)
    )
    status = get_dict_value(
        "status",
        get_dict_value("data", structured_data)
    )
    if status != "success":
        state = click.style("failed", bold=True, fg='red')
    else:
        state = click.style("passed", bold=True, fg='green')
    click.echo(test_id.label + ": " + state)
    if summary_style != "minimal":
        _print_command_line(structured_data, test_id)


def _print_command_line(test, test_id):
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
    "handle_reports_and_summary",
]
