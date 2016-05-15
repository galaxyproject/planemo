from .test_utils import (
    CliTestCase,
    skip_unless_environ,
)


class DatabaseCommandsTestCase(CliTestCase):

    @skip_unless_environ("PLANEMO_ENABLE_POSTGRES_TESTS")
    def test_database_commands(self):
        with self._isolate():
            result = self._check_exit_code(["database_list"])
            assert "test123" not in result.output
            self._check_exit_code(["database_create", "test123"])
            result = self._check_exit_code(["database_list"])
            assert "test123" in result.output
            self._check_exit_code(["database_delete", "test123"])
