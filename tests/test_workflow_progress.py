import time
from io import StringIO

from rich.console import Console

from planemo.galaxy.invocations.progress import (
    WorkflowProgress,
    WorkflowProgressDisplay,
)
from planemo.galaxy.invocations.progress_display import DisplayConfiguration

STEP_NEW = {"state": "new", "id": "1"}
STEP_SCHEDULED = {"state": "scheduled", "id": "1"}
STEP_COMPLETED = {"state": "completed", "id": "1"}
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

        workflow_progress.handle_invocation({"state": "completed", "steps": [STEP_NEW]}, {"states": {"new": 1}})
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


def test_workflow_progress_completed_invocation_state():
    """Test that 'completed' invocation state is treated like 'scheduled'."""
    completed_invocation = {"state": "completed", "steps": [STEP_COMPLETED]}

    with WorkflowProgress(DisplayConfiguration()) as workflow_progress:
        workflow_progress.handle_invocation(completed_invocation, {"states": {"ok": 2}})
        assert workflow_progress.invocation_scheduling_terminal
        assert workflow_progress.jobs_terminal
        assert workflow_progress.terminal


def test_workflow_progress_completed_step_state_counting():
    """Test that 'completed' step states are counted alongside 'scheduled' for progress."""
    mixed_steps = [STEP_SCHEDULED, STEP_COMPLETED, STEP_NEW]

    with WorkflowProgress(DisplayConfiguration()) as workflow_progress:
        workflow_progress.handle_invocation({"state": "scheduled", "steps": mixed_steps}, {"states": {"ok": 2}})
        # Both 'scheduled' and 'completed' steps should be counted
        num_scheduled = (workflow_progress.step_states.get("scheduled") or 0) + (
            workflow_progress.step_states.get("completed") or 0
        )
        assert num_scheduled == 2


def test_format_job_error_details_escapes_rich_markup():
    # Regression: failed-job stderr containing strings like Java's Arrays.toString
    # output (e.g. ImageJ --debug) starts with "[/tmp/..." which rich's markup parser
    # treats as a closing tag and rejects with MarkupError, masking the real failure.
    bracketed_argv = (
        "[/tmp/shed_dir/toolshed.g2.bx.psu.edu/repos/imgteam/imagej2_analyze_particles_binary/"
        "862af85a50ec/imagej2_analyze_particles_binary/imagej2_analyze_particles_binary_jython_script.py, "
        "/tmp/job/outputs/dataset_cc.dat, /tmp/files/dataset_99.dat, "
        "yes, 50-Infinity, 0.0, 1.0, Outlines, no, no, no, output.tiff, tiff, ]"
    )
    job_details = {
        "exit_code": 1,
        "tool_id": "toolshed.g2.bx.psu.edu/repos/imgteam/imagej2_analyze_particles_binary/"
        "imagej2_analyze_particles_binary/20240614+galaxy0",
        "command_line": "ImageJ --jython /tmp/shed_dir/script.py /tmp/input.dat",
        "stdout": "",
        "stderr": bracketed_argv,
    }

    with WorkflowProgress(DisplayConfiguration()) as workflow_progress:
        lines = workflow_progress._format_job_error_details("JOB_ID_42", job_details)

    # Render through a Console with markup parsing on (the default) to ensure
    # the assembled output does not blow up rich.
    console = Console(file=StringIO(), force_terminal=False, width=200)
    console.print("\n".join(lines))  # would raise MarkupError before the fix
    output = console.file.getvalue()

    # Content must be preserved verbatim, not stripped.
    assert "/tmp/shed_dir/" in output
    assert "output.tiff" in output
    assert "20240614+galaxy0" in output
