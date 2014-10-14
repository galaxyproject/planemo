#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_planemo
----------------------------------

Tests for `planemo` module.
"""

from click.testing import CliRunner

import unittest

from planemo import cli


# More information on testing click applications at following link.
# http://click.pocoo.org/3/testing/#basic-testing
class TestPlanemo(unittest.TestCase):

    def test_commands_have_help(self):
        commands = cli.list_cmds()
        runner = CliRunner()
        for command in commands:
            planemo_cli = cli.planemo
            result = runner.invoke(planemo_cli, [command, "--help"])
            if result.exit_code != 0:
                message = "Planemo command %s has invalid --help." % command
                raise AssertionError(message)

    def test_responds_to_desired_commands(self):
        commands = cli.list_cmds()

        def assert_responds_to(command):
            assert command in commands, "No command %s" % command

        assert_responds_to("docker_shell")
        assert_responds_to("docker_build")
        assert_responds_to("brew_init")
        assert_responds_to("brew")
        assert_responds_to("brew_env")
        assert_responds_to("config_init")
        assert_responds_to("lint")
        assert_responds_to("project_init")
        assert_responds_to("shed_upload")
        assert_responds_to("serve")
        assert_responds_to("syntax")
        assert_responds_to("test")
        assert_responds_to("travis_before_install")
        assert_responds_to("travis_init")


if __name__ == '__main__':
    unittest.main()
