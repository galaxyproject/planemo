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
