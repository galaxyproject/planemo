"""Provide abstractions over click testing of the app and unittest."""
from __future__ import print_function

import contextlib
import os
import shutil
import traceback

from sys import version_info
from tempfile import mkdtemp

from click.testing import CliRunner
from galaxy.tools.deps.commands import which

from planemo import cli
from planemo import io
from planemo import shed
from planemo.config import PLANEMO_CONFIG_ENV_PROP

from .shed_app_test_utils import (
    mock_shed,
    setup_mock_shed,
)

if version_info < (2, 7):
    from unittest2 import TestCase, skip
    PRE_PYTHON_27 = True
else:
    from unittest import TestCase, skip
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


# More information on testing click applications at following link.
# http://click.pocoo.org/3/testing/#basic-testing
class CliTestCase(TestCase):

    def setUp(self):  # noqa
        self._runner = CliRunner()
        self._home = mkdtemp()
        self._old_config = os.environ.get(PLANEMO_CONFIG_ENV_PROP, None)
        os.environ[PLANEMO_CONFIG_ENV_PROP] = self.planemo_yaml_path

    def tearDown(self):  # noqa
        if self._old_config:
            os.environ[PLANEMO_CONFIG_ENV_PROP] = self._old_config
        else:
            del os.environ[PLANEMO_CONFIG_ENV_PROP]
        shutil.rmtree(self._home)

    @property
    def planemo_yaml_path(self):
        return os.path.join(self._home, ".planemo.yml")

    @property
    def _cli(self):
        return cli

    def _isolate(self):
        return self._runner.isolated_filesystem()

    def _invoke(self, command_list):
        planemo_cli = self._cli.planemo
        return self._runner.invoke(planemo_cli, command_list)

    def _check_exit_code(self, command_list, exit_code=0):
        expected_exit_code = exit_code
        result = self._invoke(command_list)
        print("Command list output is [%s]" % result.output)
        result_exit_code = result.exit_code
        if result_exit_code != expected_exit_code:
            message = EXIT_CODE_MESSAGE % (
                " ".join(command_list),
                result_exit_code,
                expected_exit_code,
                result.output,
            )
            if result.exception:
                message += " Exception [%s], " % str(result.exception)
                exc_type, exc_value, exc_traceback = result.exc_info
                tb = traceback.format_exception(exc_type, exc_value,
                                                exc_traceback)
                message += "Traceback [%s]" % tb
            raise AssertionError(message)
        return result

    @contextlib.contextmanager
    def _isolate_repo(self, name):
        with self._isolate() as f:
            self._copy_repo(name, f)
            yield f

    def _copy_repo(self, name, dest):
        repo = os.path.join(TEST_REPOS_DIR, name)
        io.shell("cp -r '%s/.' '%s'" % (repo, dest))

    @property
    def test_context(self):
        return test_context()


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
        shutil.rmtree(self.temp_directory)


class TempDirectoryContext(object):
    def __init__(self):
        self.temp_directory = mkdtemp()

    def __enter__(self):
        return self

    def __exit__(self, type, value, tb):
        shutil.rmtree(self.temp_directory)


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


def test_context():
    context = cli.Context()
    context.planemo_directory = "/tmp/planemo-test-workspace"
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


# TODO: everything should be considered "exported".
__all__ = (
    "assert_exists",
    "TestCase",
    "CliTestCase",
)
