""" Provide abstractions over click testing of the
app and unittest.
"""
import contextlib
import os
from tempfile import mkdtemp
import shutil
from sys import version_info

from click.testing import CliRunner

from planemo import cli
from planemo import shed
from planemo import io

from galaxy.tools.deps.commands import which
from .shed_app_test_utils import (
    mock_shed,
    setup_mock_shed,
)

if version_info < (2, 7):
    from unittest2 import TestCase, skip
else:
    from unittest import TestCase, skip

TEST_DIR = os.path.dirname(__file__)
TEST_DATA_DIR = os.path.join(TEST_DIR, "data")
TEST_REPOS_DIR = os.path.join(TEST_DATA_DIR, "repos")
EXIT_CODE_MESSAGE = ("Planemo command [%s] resulted in unexpected exit code "
                     "[%s], expected exit code [%s]]. Command output [%s]")


# More information on testing click applications at following link.
# http://click.pocoo.org/3/testing/#basic-testing
class CliTestCase(TestCase):

    def setUp(self):  # noqa
        self._runner = CliRunner()

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
        result_exit_code = result.exit_code
        if result_exit_code != expected_exit_code:
            message = EXIT_CODE_MESSAGE % (
                " ".join(command_list),
                result_exit_code,
                expected_exit_code,
                result.output,
            )
            if result.exception:
                message += " Exception [%s]." % str(result.exception)
            raise AssertionError(message)
        return result

    @contextlib.contextmanager
    def _isolate_repo(self, name):
        with self._isolate() as f:
            repo = os.path.join(TEST_REPOS_DIR, name)
            io.shell("cp -r '%s'/. '%s'" % (repo, f))
            yield f


class CliShedTestCase(CliTestCase):

    def setUp(self):  # noqa
        super(CliShedTestCase, self).setUp()
        self.mock_shed = setup_mock_shed()

    def tearDown(self):  # noqa
        super(CliShedTestCase, self).tearDown()
        self.mock_shed.shutdown()

    def _shed_args(self, read_only=False):
        args = [
            "--shed_target", self.mock_shed.url,
        ]
        if not read_only:
            args.extend(["--shed_key", "ignored"])
        return args


@contextlib.contextmanager
def mock_shed_client():
    with mock_shed() as mock_shed_obj:
        yield shed.tool_shed_client(shed_target=mock_shed_obj.url)


class TempDirectoryTestCase(TestCase):

    def setUp(self):  # noqa
        self.temp_directory = mkdtemp()

    def tearDown(self):  # noqa
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


__all__ = [
    "TestCase",
    "CliTestCase",
]
