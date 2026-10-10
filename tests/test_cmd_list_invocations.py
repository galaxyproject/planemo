"""Tests for the planemo ``list_invocations`` command."""

import json
from unittest.mock import (
    call,
    Mock,
    patch,
)

import pytest

from planemo.galaxy.api import get_invocations
from .test_utils import CliTestCase


def _invocation(number, workflow_id="workflow-1"):
    return {
        "id": f"invocation-{number}",
        "workflow_id": workflow_id,
        "history_id": f"history-{number}",
    }


def test_get_invocations_limits_final_page():
    gi = Mock()
    gi.invocations.get_invocations.side_effect = [
        [_invocation(i) for i in range(20)],
        [_invocation(i) for i in range(20, 25)],
    ]
    gi.invocations.get_invocation_summary.return_value = {"states": {"ok": 1}}

    invocations = get_invocations(gi, "workflow-1", max_items=25, offset_items=5)

    assert len(invocations) == 25
    assert gi.invocations.get_invocations.call_args_list == [
        call("workflow-1", limit=20, offset=5),
        call("workflow-1", limit=5, offset=25),
    ]


def test_get_invocations_without_workflow_uses_instance_endpoint():
    gi = Mock()
    gi.invocations.get_invocations.return_value = [_invocation(1)]
    gi.invocations.get_invocation_summary.return_value = {"states": {"ok": 1}}

    assert list(get_invocations(gi, None, instance=True, max_items=1)) == ["invocation-1"]
    gi.invocations.get_invocations.assert_called_once_with(instance=True, limit=1, offset=0)


@pytest.mark.parametrize(
    "arguments",
    [
        {"max_items": -1},
        {"offset_items": -1},
        {"items_per_request": 0},
    ],
)
def test_get_invocations_rejects_invalid_pagination(arguments):
    with pytest.raises(ValueError):
        get_invocations(Mock(), None, **arguments)


class CmdListInvocationsTestCase(CliTestCase):
    def _list_invocations(self, invocations, *extra_args):
        profile = {
            "galaxy_url": "https://galaxy.example.org/",
            "galaxy_admin_key": None,
            "galaxy_user_key": "test-key",
        }
        galaxy_instance = Mock()
        galaxy_instance.workflows.show_workflow.side_effect = lambda workflow_id, **kwds: {
            "id": workflow_id,
            "name": f"Workflow {workflow_id}",
        }
        get_invocations_mock = Mock(return_value=invocations)
        runnable = Mock(is_remote_workflow_uri=True)
        with (
            patch("planemo.commands.cmd_list_invocations.profiles.ensure_profile", return_value=profile),
            patch("planemo.commands.cmd_list_invocations.gi", return_value=galaxy_instance),
            patch("planemo.commands.cmd_list_invocations.get_invocations", get_invocations_mock),
            patch("planemo.commands.cmd_list_invocations.for_runnable_identifier", return_value=runnable),
            patch("planemo.commands.cmd_list_invocations.remote_runnable_to_workflow_id", return_value="workflow-1"),
        ):
            result = self._runner.invoke(
                self._cli.planemo,
                ["list_invocations", "--profile", "test", *extra_args],
            )
        assert result.exit_code == 0, result.output
        return result, galaxy_instance, get_invocations_mock

    def test_raw_lists_invocations_without_workflow_identifier(self):
        invocations = {
            "invocation-1": {
                "states": {"ok": 1},
                "workflow_id": "workflow-1",
                "history_id": "history-1",
            }
        }

        result, galaxy_instance, get_invocations_mock = self._list_invocations(invocations, "--raw")

        assert json.loads(result.output) == invocations
        get_invocations_mock.assert_called_once_with(
            gi=galaxy_instance,
            workflow_id=None,
            instance=True,
            max_items=100,
            offset_items=0,
        )
        galaxy_instance.workflows.show_workflow.assert_not_called()

    def test_workflow_identifier_and_pagination_are_forwarded(self):
        result, galaxy_instance, get_invocations_mock = self._list_invocations(
            {},
            "workflow-alias",
            "--raw",
            "--max-items",
            "12",
            "--offset-items",
            "3",
        )

        assert result.exit_code == 0
        get_invocations_mock.assert_called_once_with(
            gi=galaxy_instance,
            workflow_id="workflow-1",
            instance=True,
            max_items=12,
            offset_items=3,
        )

    def test_table_groups_invocations_by_workflow(self):
        invocations = {
            "invocation-1": {
                "states": {"ok": 1, "custom_state": 2},
                "workflow_id": "workflow-1",
                "history_id": "history-1",
            },
            "invocation-2": {
                "states": {"running": 1},
                "workflow_id": "workflow-2",
                "history_id": "history-2",
            },
        }

        result, galaxy_instance, _ = self._list_invocations(invocations)

        assert "Looking for invocations for all workflows" in result.output
        assert "Workflow workflow-1" in result.output
        assert "Workflow workflow-2" in result.output
        assert "2 jobs custom_state" in result.output
        assert "2 invocations found." in result.output
        assert galaxy_instance.workflows.show_workflow.call_count == 2
