import time

from planemo.galaxy.invocations.progress import (
    WorkflowProgress,
    WorkflowProgressDisplay,
)
from planemo.galaxy.invocations.progress_display import DisplayConfiguration

STEP_NEW = {"state": "new", "id": "1"}
STEP_SCHEDULED = {"state": "scheduled", "id": "1"}
SLEEP = 0.8


def test_workflow_progress_typical():
    with WorkflowProgressDisplay("myid12345abcde") as live:
        new_steps = [STEP_NEW, STEP_NEW, STEP_NEW, STEP_NEW]
        one_scheduled_steps = [STEP_SCHEDULED, STEP_NEW, STEP_NEW, STEP_NEW]
        two_scheduled_steps = [STEP_SCHEDULED, STEP_SCHEDULED, STEP_NEW, STEP_NEW]
        three_scheduled_steps = [STEP_SCHEDULED, STEP_SCHEDULED, STEP_SCHEDULED, STEP_NEW]
        all_scheduled_steps = [STEP_SCHEDULED, STEP_SCHEDULED, STEP_SCHEDULED, STEP_SCHEDULED]

        state_pairs = [
            ({"state": "new"}, {}, None),
            ({"state": "ready", "steps": new_steps}, {}, None),
            ({"state": "ready", "steps": one_scheduled_steps}, {"states": {"new": 1}}, None),
            ({"state": "ready", "steps": two_scheduled_steps}, {"states": {"new": 2}}, None),
            ({"state": "ready", "steps": two_scheduled_steps}, {"states": {"new": 1, "running": 1}}, None),
            ({"state": "ready", "steps": two_scheduled_steps}, {"states": {"new": 1, "ok": 1}}, None),
            ({"state": "ready", "steps": two_scheduled_steps}, {"states": {"ok": 2}}, None),
            ({"state": "scheduled", "steps": three_scheduled_steps}, {"states": {"ok": 2, "new": 3}}, None),
            (
                {"state": "scheduled", "steps": three_scheduled_steps},
                {"states": {"ok": 2, "running": 1, "new": 2}},
                None,
            ),
            (
                {"state": "scheduled", "steps": three_scheduled_steps},
                {"states": {"ok": 3, "running": 1, "new": 1}},
                None,
            ),
            ({"state": "scheduled", "steps": three_scheduled_steps}, {"states": {"ok": 4, "running": 1}}, None),
            ({"state": "scheduled", "steps": three_scheduled_steps}, {"states": {"ok": 5}}, None),
            ({"state": "scheduled", "steps": all_scheduled_steps}, {"states": {"ok": 5}}, None),
            ({"state": "scheduled", "steps": all_scheduled_steps}, {"states": {"ok": 5}}, None),
            ({"state": "ready", "steps": new_steps}, {}, "abcde123456"),
            ({"state": "ready", "steps": one_scheduled_steps}, {"states": {"new": 1}}, "abcde123456"),
            ({"state": "ready", "steps": two_scheduled_steps}, {"states": {"new": 2}}, "abcde123456"),
            ({"state": "ready", "steps": two_scheduled_steps}, {"states": {"new": 1, "running": 1}}, "abcde123456"),
            ({"state": "ready", "steps": two_scheduled_steps}, {"states": {"new": 1, "ok": 1}}, "abcde123456"),
            ({"state": "ready", "steps": two_scheduled_steps}, {"states": {"ok": 2}}, "abcde123456"),
            ({"state": "scheduled", "steps": three_scheduled_steps}, {"states": {"ok": 2, "new": 3}}, "abcde123456"),
            (
                {"state": "scheduled", "steps": three_scheduled_steps},
                {"states": {"ok": 2, "running": 1, "new": 2}},
                "abcde123456",
            ),
            (
                {"state": "scheduled", "steps": three_scheduled_steps},
                {"states": {"ok": 3, "running": 1, "new": 1}},
                "abcde123456",
            ),
            (
                {"state": "scheduled", "steps": three_scheduled_steps},
                {"states": {"ok": 4, "running": 1}},
                "abcde123456",
            ),
            ({"state": "scheduled", "steps": three_scheduled_steps}, {"states": {"ok": 5}}, "abcde123456"),
            ({"state": "scheduled", "steps": all_scheduled_steps}, {"states": {"ok": 5}}, "abcde123456"),
            ({"state": "scheduled", "steps": all_scheduled_steps}, {"states": {"ok": 5}}, "abcde123456"),
        ]
        i = 0
        for invocation, job_states_summary, subworkflow_invocation_id in state_pairs:
            i = i + 1
            if subworkflow_invocation_id is None:
                live.handle_invocation(invocation, job_states_summary)
            else:
                invocation = invocation.copy()
                invocation["id"] = subworkflow_invocation_id
                live.handle_subworkflow_invocation(invocation, job_states_summary)
            time.sleep(SLEEP)


def test_workflow_progress_scheduling_state_handling():
    with WorkflowProgress(DisplayConfiguration()) as workflow_progress:
        workflow_progress.handle_invocation({"state": "new", "steps": [STEP_NEW]}, {"states": {"new": 1}})
        assert not workflow_progress.invocation_scheduling_terminal

        workflow_progress.handle_invocation({"state": "ready", "steps": [STEP_NEW]}, {"states": {"new": 1}})
        assert not workflow_progress.invocation_scheduling_terminal

        workflow_progress.handle_invocation({"state": "cancelling", "steps": [STEP_NEW]}, {"states": {"new": 1}})
        assert not workflow_progress.invocation_scheduling_terminal

        workflow_progress.handle_invocation(
            {"state": "requires_materialization", "steps": [STEP_NEW]}, {"states": {"new": 1}}
        )
        assert not workflow_progress.invocation_scheduling_terminal

        workflow_progress.handle_invocation({"state": "scheduled", "steps": [STEP_NEW]}, {"states": {"new": 1}})
        assert workflow_progress.invocation_scheduling_terminal

        workflow_progress.handle_invocation({"state": "cancelled", "steps": [STEP_NEW]}, {"states": {"new": 1}})
        assert workflow_progress.invocation_scheduling_terminal

        workflow_progress.handle_invocation({"state": "failed", "steps": [STEP_NEW]}, {"states": {"new": 1}})
        assert workflow_progress.invocation_scheduling_terminal


# From Galaxy:
# NEW = "new"
# RESUBMITTED = "resubmitted"
# UPLOAD = "upload"
# WAITING = "waiting"
# QUEUED = "queued"
# RUNNING = "running"
# OK = "ok"
# ERROR = "error"
# FAILED = "failed"
# PAUSED = "paused"
# DELETING = "deleting"
# DELETED = "deleted"
# STOPPING = "stop"
# STOPPED = "stopped"
# SKIPPED = "skipped"
def test_workflow_progress_job_state_handling():
    scheduled_invocation = {"state": "scheduled", "steps": [STEP_SCHEDULED]}

    with WorkflowProgress(DisplayConfiguration()) as workflow_progress:
        workflow_progress.handle_invocation(scheduled_invocation, {"states": {"new": 1}})
        assert not workflow_progress.jobs_terminal

        workflow_progress.handle_invocation(scheduled_invocation, {"states": {"new": 1, "ok": 2}})
        assert not workflow_progress.jobs_terminal

        workflow_progress.handle_invocation(scheduled_invocation, {"states": {"ok": 2}})
        assert workflow_progress.jobs_terminal

        workflow_progress.handle_invocation(scheduled_invocation, {"states": {"ok": 2, "paused": 1}})
        assert workflow_progress.jobs_terminal

        workflow_progress.handle_invocation(scheduled_invocation, {"states": {"ok": 2, "paused": 1, "new": 1}})
        assert not workflow_progress.jobs_terminal
