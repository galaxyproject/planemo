from planemo.database.postgres_docker import stop_postgres_docker
from .test_utils import (
    CliTestCase,
    skip_unless_environ,
    skip_unless_executable,
)


class ProfileCommandsTestCase(CliTestCase):
    def _profile_commands(self, database_type="postgres"):
        with self._isolate():
            result = self._check_exit_code(["profile_list"])
            assert "profile1234" not in result.output
            self._check_exit_code(["profile_create", "profile1234", "--database_type", database_type])
            result = self._check_exit_code(["profile_list"])
            assert "profile1234" in result.output
            self._check_exit_code(["profile_delete", "profile1234"])
            result = self._check_exit_code(["profile_list"])
            assert "profile1234" not in result.output, result.output

    @skip_unless_environ("PLANEMO_ENABLE_POSTGRES_TESTS")
    @skip_unless_executable("psql")
    def test_profile_commands(self):
        self._profile_commands()

    @skip_unless_executable("docker")
    def test_profile_commands_docker(self):
        try:
            self._profile_commands(database_type="postgres_docker")
        finally:
            stop_postgres_docker()
