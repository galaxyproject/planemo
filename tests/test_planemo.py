#!/usr/bin/env python

"""
test_planemo
----------------------------------

Tests for `planemo` module.
"""

from .test_utils import CliTestCase


class TestPlanemo(CliTestCase):
    def test_commands_have_help(self):
        commands = self._cli.list_cmds()
        for command in commands:
            self._check_exit_code([command, "--help"])

    def test_responds_to_desired_commands(self):
        commands = self._cli.list_cmds()

        def assert_responds_to(command):
            assert command in commands, "No command %s" % command

        assert_responds_to("config_init")
        assert_responds_to("create_gist")
        assert_responds_to("docker_build")
        assert_responds_to("docker_shell")
        assert_responds_to("lint")
        assert_responds_to("normalize")
        assert_responds_to("project_init")
        assert_responds_to("serve")
        assert_responds_to("shed_diff")
        assert_responds_to("shed_download")
        assert_responds_to("shed_upload")
        assert_responds_to("syntax")
        assert_responds_to("test")
        assert_responds_to("tool_init")

    def test_planemo_version_command(self):
        self._check_exit_code(["--version"])

    def test_planemo_help_command(self):
        self._check_exit_code(["--help"])

    def test_unknown_command_error(self):
        """Test that unknown commands produce helpful error messages."""
        result = self._check_exit_code(["profile_init", "--help"], exit_code=self.non_zero_exit_code)
        # Should show error message
        assert "No such command 'profile_init'" in result.output
        # Should suggest similar commands
        assert "Did you mean one of these?" in result.output
        assert "profile_create" in result.output

    def test_unknown_command_with_suggestions(self):
        """Test that typos in commands suggest similar alternatives."""
        result = self._check_exit_code(["workflow_tst", "--help"], exit_code=self.non_zero_exit_code)
        assert "No such command 'workflow_tst'" in result.output
        # Should suggest workflow commands since they start with "workflow_"
        assert "Did you mean one of these?" in result.output

    def test_completely_unknown_command(self):
        """Test that completely invalid commands show basic error."""
        result = self._check_exit_code(["foobar", "--help"], exit_code=self.non_zero_exit_code)
        assert "No such command 'foobar'" in result.output
        assert "planemo --help" in result.output
