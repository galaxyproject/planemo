#!/usr/bin/env python

import json

from planemo.cli import INTERNAL_COMMANDS

from .test_utils import CliTestCase


class TestCliMetadata(CliTestCase):
    def test_command_names_match_modules(self):
        for command_name in self._cli.list_cmds():
            command = self._cli.name_to_command(command_name)
            assert command.name == command_name

    def test_internal_command_policy(self):
        assert INTERNAL_COMMANDS == ["create_gist", "shed_download"]

    def test_cli_metadata_outputs_json(self):
        result = self._check_exit_code(["cli_metadata", "--format", "json"])
        metadata = json.loads(result.output)

        assert metadata["schema_version"] == "0.1"
        assert metadata["program"] == "planemo"
        assert metadata["aliases"] == {"l": "lint", "o": "open", "s": "serve", "t": "test"}

        command_names = {command["name"] for command in metadata["commands"]}
        assert "test" in command_names
        assert "run" in command_names
        assert "workflow_test_init" in command_names
        assert "workflow_test_on_invocation" in command_names
        assert "invocation_download" in command_names
        assert "create_gist" not in command_names
        assert "shed_download" not in command_names

    def test_cli_metadata_includes_internal_when_requested(self):
        result = self._check_exit_code(["cli_metadata", "--format", "json", "--include-internal"])
        metadata = json.loads(result.output)
        command_names = {command["name"] for command in metadata["commands"]}

        assert "create_gist" in command_names
        assert "shed_download" in command_names

    def test_cli_metadata_for_test_command(self):
        result = self._check_exit_code(["cli_metadata", "--format", "json", "--command", "test"])
        metadata = json.loads(result.output)

        assert metadata["name"] == "test"
        params = {param["name"]: param for param in metadata["params"]}
        assert "failed" in params
        assert "test_index" in params
        assert params["test_output_json"]["default"] == "tool_test_output.json"
        assert "--test_output_json" in params["test_output_json"]["opts"]
        assert "cwltool" in params["engine"]["type"]["choices"]

    def test_cli_metadata_for_run_command(self):
        result = self._check_exit_code(["cli_metadata", "--format", "json", "--command", "run"])
        metadata = json.loads(result.output)

        params = {param["name"]: param for param in metadata["params"]}
        assert "runnable_identifier" in params
        assert "job_path" in params
        assert "export_format" in params
        assert "engine" in params
        assert "output_json" in params
        assert "--output_json" in params["output_json"]["opts"]

    def test_cli_metadata_for_invocation_download_command(self):
        result = self._check_exit_code(["cli_metadata", "--format", "json", "--command", "invocation_download"])
        metadata = json.loads(result.output)

        params = {param["name"]: param for param in metadata["params"]}
        assert "invocation_id" in params
        assert "output_json" in params
        assert "output_json_path_type" in params
        assert params["output_json_path_type"]["type"]["choices"] == ["relative", "absolute"]
        assert "--no_ignore_missing_output" in params["ignore_missing_output"]["secondary_opts"]

    def test_cli_metadata_for_workflow_test_on_invocation_command(self):
        result = self._check_exit_code(["cli_metadata", "--format", "json", "--command", "workflow_test_on_invocation"])
        metadata = json.loads(result.output)

        arguments = [param for param in metadata["params"] if param["kind"] == "argument"]
        assert [argument["name"] for argument in arguments] == ["path", "invocation_id"]
        assert [argument["human_readable_name"] for argument in arguments] == ["TEST.YML", "INVOCATION_ID"]
