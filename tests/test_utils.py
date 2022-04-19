"""Provide abstractions over click testing of the app and unittest."""
from __future__ import print_function

import contextlib
import functools
import os
import re
import shutil
import signal
import traceback
from concurrent.futures import as_completed, ThreadPoolExecutor
from sys import version_info
from tempfile import mkdtemp
from unittest import skip, TestCase

import psutil
import py.code
import pytest
from click.testing import CliRunner
from galaxy.util import asbool, unicodify, which

from planemo import cli
from planemo import io
from planemo import shed
from planemo.config import PLANEMO_CONFIG_ENV_PROP
from planemo.galaxy.ephemeris_sleep import sleep, SleepCondition
from .shed_app_test_utils import (
    mock_shed,
    setup_mock_shed,
)

PRE_PYTHON_27 = False
if version_info[0] == 2 and version_info[1] >= 7:
    PYTHON_27 = True
else:
    PYTHON_27 = False

TEST_DIR = os.path.dirname(__file__)
TEST_DATA_DIR = os.path.join(TEST_DIR, "data")
TEST_REPOS_DIR = os.path.join(TEST_DATA_DIR, "repos")
TEST_RECIPES_DIR = os.path.join(TEST_DATA_DIR, "recipes")
TEST_TOOLS_DIR = os.path.join(TEST_DATA_DIR, "tools")
PROJECT_TEMPLATES_DIR = os.path.join(TEST_DIR, os.path.pardir, "project_templates")
EXIT_CODE_MESSAGE = ("Planemo command [%s] resulted in unexpected exit code "
                     "[%s], expected exit code [%s]]. Command output [%s]")
CWL_DRAFT3_DIR = os.path.join(PROJECT_TEMPLATES_DIR, "cwl_draft3_spec")
NON_ZERO_EXIT_CODE = object()


class MarkGenerator(object):

    def __getattr__(self, name):
        return getattr(pytest.mark, name)


mark = MarkGenerator()


# More information on testing click applications at following link.
# http://click.pocoo.org/3/testing/#basic-testing
class CliTestCase(TestCase):
    non_zero_exit_code = NON_ZERO_EXIT_CODE

    def setUp(self):  # noqa
        self._runner = CliRunner()
        self._home = mkdtemp()
        self._old_config = os.environ.get(PLANEMO_CONFIG_ENV_PROP, None)
        self._futures = []
        self._port = None
        os.environ[PLANEMO_CONFIG_ENV_PROP] = self.planemo_yaml_path

    def tearDown(self):  # noqa
        for future in self._futures:
            future.cancel()
            try:
                future.result(timeout=10)
            except Exception as e:
                print("Failed to dispose of future driven resource [%s]" % e)
        if self._port:
            kill_process_on_port(self._port)
        if self._old_config:
            os.environ[PLANEMO_CONFIG_ENV_PROP] = self._old_config
        else:
            del os.environ[PLANEMO_CONFIG_ENV_PROP]
        safe_rmtree(self._home)

    @property
    def planemo_yaml_path(self):
        return os.path.join(self._home, ".planemo.yml")

    @property
    def _cli(self):
        return cli

    def _isolate(self):
        return self._runner.isolated_filesystem()

    def _check_exit_code(self, command_list, exit_code=0):
        return check_exit_code(self._runner, command_list, exit_code=exit_code)

    @contextlib.contextmanager
    def _isolate_repo(self, name):
        with self._isolate() as f:
            self._copy_repo(name, f)
            yield f

    @contextlib.contextmanager
    def _isolate_with_test_data(self, relative_path):
        with self._isolate() as f:
            repo = os.path.join(TEST_DATA_DIR, relative_path)
            self._copy_directory(repo, f)
            yield f

    def _copy_repo(self, name, dest):
        repo = os.path.join(TEST_REPOS_DIR, name)
        self._copy_directory(repo, dest)

    def _copy_directory(self, path, dest):
        io.shell(['cp', '-r', "%s/." % path, dest])

    @property
    def test_context(self):
        return create_test_context()


class CliShedTestCase(CliTestCase):

    def setUp(self):  # noqa
        super(CliShedTestCase, self).setUp()
        self.mock_shed = setup_mock_shed()

    def tearDown(self):  # noqa
        super(CliShedTestCase, self).tearDown()
        self.mock_shed.shutdown()

    def _shed_create(self, recursive=False):
        create_command = ["shed_create"]
        if recursive:
            create_command.append("-r")
        create_command.extend(self._shed_args())
        self._check_exit_code(create_command)

    def _shed_args(self, read_only=False):
        args = [
            "--shed_target", self.mock_shed.url,
        ]
        if not read_only:
            args.extend(["--shed_key", "ignored"])
        return args

    def _print_shed_info(self):
        print(self._repositories)

    @property
    def _tsi(self):
        tsi = shed.tool_shed_client(
            None,
            shed_target=self.mock_shed.url,
            key="ignored",
        )
        return tsi

    @property
    def _repositories(self):
        return self._tsi.repositories.get_repositories()

    def repository_by_name(self, name):
        return [r for r in self._repositories if r["name"] == name][0]


@contextlib.contextmanager
def mock_shed_context():
    with mock_shed() as mock_shed_obj:
        yield shed.get_shed_context(shed_target=mock_shed_obj.url)


class TempDirectoryTestCase(TestCase):

    def setUp(self):  # noqa
        self.temp_directory = mkdtemp()

    def tearDown(self):  # noqa
        safe_rmtree(self.temp_directory)


class TempDirectoryContext(object):
    def __init__(self):
        self.temp_directory = mkdtemp()

    def __enter__(self):
        return self

    def __exit__(self, type, value, tb):
        safe_rmtree(self.temp_directory)


def skip_unless_environ(var):
    if var in os.environ:
        return lambda func: func
    template = "Environment variable %s not found, dependent test skipped."
    return skip(template % var)


def skip_if_environ(var):
    if var not in os.environ:
        return lambda func: func
    template = "Environment variable %s set, dependent test skipped."
    return skip(template % var)


def skip_unless_module(module):
    available = True
    try:
        __import__(module)
    except ImportError:
        available = False
    if available:
        return lambda func: func
    template = "Module %s could not be loaded, dependent test skipped."
    return skip(template % module)


def skip_unless_executable(executable):
    if which(executable):
        return lambda func: func
    return skip("PATH doesn't contain executable %s" % executable)


def skip_unless_python_2_7():
    if PYTHON_27:
        return lambda func: func
    return skip("Python 2.7 required for test.")


def target_galaxy_branch():
    return os.environ.get("PLANEMO_TEST_GALAXY_BRANCH", "master")


# Taken from Galaxy's test/unit/tools/test_tool_deps.py
@contextlib.contextmanager
def modify_environ(values, remove=[]):
    """
    Modify the environment for a test, adding/updating values in dict `values` and
    removing any environment variables mentioned in list `remove`.
    """
    new_keys = set(values.keys()) - set(os.environ.keys())
    old_environ = os.environ.copy()
    try:
        os.environ.update(values)
        for to_remove in remove:
            try:
                del os.environ[remove]
            except KeyError:
                pass
        yield
    finally:
        os.environ.update(old_environ)
        for key in new_keys:
            del os.environ[key]


def run_verbosely():
    return asbool(os.environ.get("PLANEMO_TEST_VERBOSE", "false"))


def create_test_context():
    context = cli.PlanemoCliContext()
    context.planemo_directory = "/tmp/planemo-test-workspace"
    context.verbose = run_verbosely()
    return context


def assert_equal(a, b):
    """Assert two things are equal."""
    assert a == b, "%s != %s" % (a, b)


def assert_exists(path):
    """Assert supplied ``path`` exists.

    Produces an informative :class:`AssertionError` if it is does not.
    """
    dir_path = os.path.dirname(path)
    msg = None
    if not os.path.exists(dir_path):
        template = "Expected path [%s] to exist, but parent absent."
        msg = template % path
    if not os.path.exists(path):
        contents = os.listdir(dir_path)
        template = "Expected path [%s] to exist. Directory contents %s."
        msg = template % (path, contents)
    if msg is not None:
        raise AssertionError(msg)


def check_exit_code(runner, command_list, exit_code=0):
    expected_exit_code = exit_code
    planemo_cli = cli.planemo
    if run_verbosely():
        print("Invoking command [%s]" % command_list)
    result = runner.invoke(planemo_cli, command_list)
    if run_verbosely():
        print("Command list output is [%s]" % result.output)
    result_exit_code = result.exit_code
    if expected_exit_code is NON_ZERO_EXIT_CODE:
        matches_expectation = result_exit_code != 0
    else:
        matches_expectation = result_exit_code == expected_exit_code
    if not matches_expectation:
        message = EXIT_CODE_MESSAGE % (
            " ".join(command_list),
            result_exit_code,
            expected_exit_code,
            result.output,
        )
        if result.exception:
            message += " Exception [%s], " % unicodify(result.exception)
            exc_type, exc_value, exc_traceback = result.exc_info
            tb = traceback.format_exception(exc_type, exc_value,
                                            exc_traceback)
            message += "Traceback [%s]" % tb
        if run_verbosely():
            print("Raising assertion error for unexpected exit code [%s]" % message)
        raise AssertionError(message)
    return result


def kill_process_on_port(port):
    # based on https://stackoverflow.com/a/20691431
    processes = []
    for proc in psutil.process_iter():
        try:
            for conns in proc.connections(kind='inet'):
                if conns.laddr.port == port:
                    proc.send_signal(signal.SIGINT)
                    processes.append(proc)
                    continue
        except Exception:
            pass


@contextlib.contextmanager
def cli_daemon_galaxy(runner, pid_file, port, command_list, exit_code=0):
    future = launch_and_wait_for_galaxy(port, check_exit_code, args=[runner, command_list, exit_code])
    yield
    io.kill_pid_file(pid_file)
    _wait_on_future_suppress_exception(future)


def launch_and_wait_for_galaxy(port, func, args=[], timeout=600, timeout_multiplier=1):
    """Run func(args) in a thread and wait on port for service.

    Service should remain up so check network a few times, this prevents
    the code that finds a free port from causing a false positive when
    detecting that the port is bound to.
    """
    target = functools.partial(func, *args)

    wait_sleep_condition = SleepCondition()

    def wait():
        effective_timeout = timeout * timeout_multiplier
        if not sleep("http://localhost:%d" % port, verbose=True, timeout=effective_timeout, sleep_condition=wait_sleep_condition):
            raise Exception('Galaxy failed to start on port %d' % port)

    executor = ThreadPoolExecutor(max_workers=2)
    try:
        target_future = executor.submit(target)
        wait_future = executor.submit(wait)
        for first in as_completed([target_future, wait_future]):
            break

        if target_future.running():
            # If wait timed out re-throw.
            wait_future.result()
            return target_future
        else:
            if target_future.exception() is not None:
                wait_future.cancel()
                wait_sleep_condition.cancel()
                _wait_on_future_suppress_exception(wait_future)
                # If the target threw an exception, rethrow.
                target_future.result()
            else:
                # Otherwise, daemon started properly. Just wait on wait.
                wait_future.result()
                return target_future
    finally:
        executor.shutdown(wait=False)


def _wait_on_future_suppress_exception(future):
    try:
        future.result(timeout=30)
    except Exception as e:
        print("Problem waiting on future %s" % e)


# From pytest-raisesregexp
class assert_raises_regexp(object):
    def __init__(self, expected_exception, regexp, *args, **kwargs):
        __tracebackhide__ = True
        self.exception = expected_exception
        self.regexp = regexp
        self.excinfo = None

        if args:
            with self:
                args[0](*args[1:], **kwargs)

    def __enter__(self):
        self.excinfo = object.__new__(py.code.ExceptionInfo)
        return self.excinfo

    def __exit__(self, exc_type, exc_val, exc_tb):
        __tracebackhide__ = True

        if exc_type is None:
            pytest.fail('DID NOT RAISE {0}'.format(self.exception))

        self.excinfo.__init__((exc_type, exc_val, exc_tb))

        if not issubclass(exc_type, self.exception):
            pytest.fail('{0} RAISED instead of {1}\n{2!r}'
                        .format(exc_type, self.exception, exc_val))

        if not re.search(self.regexp, str(exc_val)):
            pytest.fail('Pattern "{0}" not found in "{1!s}"'
                        .format(self.regexp, exc_val))

        return True


def safe_rmtree(path):
    try:
        shutil.rmtree(path)
    except Exception as e:
        print("Failed to cleanup test directory [%s]: [%s]" % (path, e))


# TODO: everything should be considered "exported".
__all__ = (
    "assert_exists",
    "cli_daemon_galaxy",
    "launch_and_wait_for_galaxy",
    "TestCase",
    "CliTestCase",
)
