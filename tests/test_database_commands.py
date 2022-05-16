from planemo.database.postgres_docker import stop_postgres_docker
from .test_utils import (
    CliTestCase,
    skip_unless_environ,
    skip_unless_executable,
)


class DatabaseCommandsTestCase(CliTestCase):
    @skip_unless_environ("PLANEMO_ENABLE_POSTGRES_TESTS")
    def _database_commands(self, database_type="postgres"):
        with self._isolate():
            result = self._check_exit_code(["database_list", "--database_type", database_type])
            assert "test1234" not in result.output
            self._check_exit_code(["database_create", "test1234", "--database_type", database_type])
            result = self._check_exit_code(["database_list", "--database_type", database_type])
            assert "test1234" in result.output
            self._check_exit_code(["database_delete", "test1234", "--database_type", database_type])
            result = self._check_exit_code(["database_list", "--database_type", database_type])
            assert "test1234" not in result.output

    @skip_unless_environ("PLANEMO_ENABLE_POSTGRES_TESTS")
    @skip_unless_executable("psql")
    def test_database_commands(self):
        self._database_commands()

    @skip_unless_environ("PLANEMO_ENABLE_POSTGRES_TESTS")
    @skip_unless_executable("docker")
    def test_database_commands_docker(self):
        try:
            self._database_commands(database_type="postgres_docker")
        finally:
            stop_postgres_docker()
