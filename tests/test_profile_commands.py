from .test_utils import (
    CliTestCase,
)


class ProfileCommandsTestCase(CliTestCase):

    def test_profile_commands(self):
        with self._isolate():
            result = self._check_exit_code(["profile_list"])
            assert "profile1234" not in result.output
            self._check_exit_code(["profile_create", "profile1234"])
            result = self._check_exit_code(["profile_list"])
            assert "profile1234" in result.output
            self._check_exit_code(["profile_delete", "profile1234"])
            result = self._check_exit_code(["profile_list"])
            assert "profile1234" not in result.output, result.output
