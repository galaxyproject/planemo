"""Provide abstractions over click testing of the app and unittest."""

import contextlib
import functools
import os
import shutil
import signal
import traceback
from concurrent.futures import (
    as_completed,
    ThreadPoolExecutor,
)
from sys import version_info
from tempfile import mkdtemp
from typing import (
    Callable,
    List,
)
from unittest import (
    skip,
    TestCase,
)

import psutil
import pytest
from click.testing import CliRunner
from galaxy.util import (
    asbool,
    unicodify,
    which,
)

from planemo import (
    cli,
    io,
    shed,
)
from planemo.config import PLANEMO_CONFIG_ENV_PROP
from planemo.galaxy.ephemeris_sleep import (
    sleep,
    SleepCondition,
)
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
TEST_TOOLS_DIR = os.path.join(TEST_DATA_DIR, "tools")
PROJECT_TEMPLATES_DIR = os.path.join(TEST_DIR, os.path.pardir, "project_templates")
CWL_DRAFT3_DIR = os.path.join(PROJECT_TEMPLATES_DIR, "cwl_draft3_spec")
NON_ZERO_EXIT_CODE = object()


class MarkGenerator:
    def __getattr__(self, name):
        return getattr(pytest.mark, name)


mark = MarkGenerator()


# More information on testing click applications at following link.
# http://click.pocoo.org/3/testing/#basic-testing
class CliTestCase(TestCase):
    non_zero_exit_code = NON_ZERO_EXIT_CODE

    def setUp(self) -> None:
        self._runner = CliRunner()
        self._home = mkdtemp()
        self._old_config = os.environ.get(PLANEMO_CONFIG_ENV_PROP, None)
        self._cleanup_hooks: List[Callable] = []
        self._port = None
        os.environ[PLANEMO_CONFIG_ENV_PROP] = self.planemo_yaml_path

    def tearDown(self):  # noqa
        for cleanup_hook in self._cleanup_hooks:
            try:
                cleanup_hook()
            except Exception as e:
                print(f"Failed to run cleanup hook: [{e}]")
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
        io.shell(["cp", "-r", f"{path}/.", dest])

    @property
    def test_context(self):
        return create_test_context()


class CliShedTestCase(CliTestCase):
    def setUp(self):  # noqa
        super().setUp()
        self.mock_shed = setup_mock_shed()

    def tearDown(self):  # noqa
        super().tearDown()
        self.mock_shed.shutdown()

    def _shed_create(self, recursive=False):
        create_command = ["shed_create"]
        if recursive:
            create_command.append("-r")
        create_command.extend(self._shed_args())
        self._check_exit_code(create_command)

    def _shed_args(self, read_only=False):
        args = [
            "--shed_target",
            self.mock_shed.url,
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


class TempDirectoryContext:
    def __init__(self):
        self.temp_directory = mkdtemp()

    def __enter__(self):
        return self

    def __exit__(self, type, value, tb):
        safe_rmtree(self.temp_directory)


def skip_unless_environ(var):
    if var in os.environ:
        return lambda func: func
    return skip(f"Environment variable {var} not found, dependent test skipped.")


def skip_if_environ(var):
    if var not in os.environ:
        return lambda func: func
    return skip(f"Environment variable {var} set, dependent test skipped.")


def skip_unless_module(module):
    available = True
    try:
        __import__(module)
    except ImportError:
        available = False
    if available:
        return lambda func: func
    return skip(f"Module {module} could not be loaded, dependent test skipped.")


def skip_unless_executable(executable):
    if which(executable):
        return lambda func: func
    return skip(f"PATH doesn't contain executable {executable}")


def target_galaxy_branch():
    return os.environ.get("PLANEMO_TEST_GALAXY_BRANCH", "master")


# Taken from Galaxy's test/unit/tool_util/util.py
@contextlib.contextmanager
def modify_environ(values, keys_to_remove=None):
    """
    Modify the environment for a test, adding/updating values in dict `values` and
    removing any environment variables mentioned in list `remove`.
    """
    old_environ = os.environ.copy()
    try:
        if values:
            os.environ.update(values)
        if keys_to_remove:
            for key in keys_to_remove:
                if key in os.environ:
                    del os.environ[key]
        yield
    finally:
        os.environ.clear()
        os.environ.update(old_environ)


def run_verbosely():
    return asbool(os.environ.get("PLANEMO_TEST_VERBOSE", "false"))


def create_test_context():
    context = cli.PlanemoCliContext()
    context.planemo_directory = "/tmp/planemo-test-workspace"
    context.verbose = run_verbosely()
    return context


def assert_equal(a, b):
    """Assert two things are equal."""
    assert a == b, f"{a} != {b}"


def assert_exists(path):
    """Assert supplied ``path`` exists.

    Produces an informative :class:`AssertionError` if it is does not.
    """
    dir_path = os.path.dirname(path)
    msg = None
    if not os.path.exists(dir_path):
        msg = f"Expected path [{path}] to exist, but parent absent."
    if not os.path.exists(path):
        contents = os.listdir(dir_path)
        msg = f"Expected path [{path}] to exist. Directory contents {contents}."
    if msg is not None:
        raise AssertionError(msg)


def check_exit_code(runner, command_list, exit_code=0):
    expected_exit_code = exit_code
    planemo_cli = cli.planemo
    if run_verbosely():
        print(f"Invoking command [{command_list}]")
    result = runner.invoke(planemo_cli, command_list)
    if run_verbosely():
        print(f"Command list output is [{result.output}]")
    result_exit_code = result.exit_code
    if expected_exit_code is NON_ZERO_EXIT_CODE:
        matches_expectation = result_exit_code != 0
    else:
        matches_expectation = result_exit_code == expected_exit_code
    if not matches_expectation:
        message = (
            f"Planemo command [{' '.join(command_list)}] resulted in unexpected exit code [{result_exit_code}], "
            f"expected exit code [{expected_exit_code}]]. Command output [{result.output}]"
        )
        if result.exception:
            message += f" Exception [{unicodify(result.exception)}], "
            exc_type, exc_value, exc_traceback = result.exc_info
            tb = traceback.format_exception(exc_type, exc_value, exc_traceback)
            message += f"Traceback [{tb}]"
        if run_verbosely():
            print(f"Raising assertion error for unexpected exit code [{message}]")
        raise AssertionError(message)
    return result


def kill_process_on_port(port):
    # based on https://stackoverflow.com/a/20691431
    processes = []
    for proc in psutil.process_iter():
        try:
            for conns in proc.connections(kind="inet"):
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


def launch_and_wait_for_galaxy(port, func, args=[], timeout=600, timeout_multiplier=1, run_as_subprocess=False):
    """Run func(args) in a thread and wait on port for service.

    Service should remain up so check network a few times, this prevents
    the code that finds a free port from causing a false positive when
    detecting that the port is bound to.
    """
    target = None
    if run_as_subprocess:
        # func is responsible for starting subprocess registering a subprocess.Popen instance for cleanup.
        # Should return immediately.
        func(*args)
    else:
        target = functools.partial(func, *args)

    wait_sleep_condition = SleepCondition()

    def wait():
        effective_timeout = timeout * timeout_multiplier
        if not sleep(
            f"http://localhost:{port}", verbose=True, timeout=effective_timeout, sleep_condition=wait_sleep_condition
        ):
            raise Exception(f"Galaxy failed to start on port {port}")

    executor = ThreadPoolExecutor(max_workers=2)
    try:
        futures = []
        if target:
            target_future = executor.submit(target)
            futures.append(target_future)
        else:
            target_future = None
        wait_future = executor.submit(wait)
        futures.append(wait_future)

        for _ in as_completed(futures):
            break

        if target_future and target_future.running():
            # If wait timed out re-throw.
            wait_future.result()
            return target_future
        else:
            if target_future and target_future.exception() is not None:
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
        print(f"Problem waiting on future {e}")


def safe_rmtree(path):
    try:
        shutil.rmtree(path)
    except Exception as e:
        print(f"Failed to cleanup test directory [{path}]: [{e}]")


# TODO: everything should be considered "exported".
__all__ = (
    "assert_exists",
    "cli_daemon_galaxy",
    "launch_and_wait_for_galaxy",
    "TestCase",
    "CliTestCase",
)
