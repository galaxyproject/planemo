import time
from typing import (
    List,
    Optional,
    Protocol,
)

from .api import (
    invocation_state_terminal,
    InvocationApi,
    InvocationJobsSummary,
    JOB_ERROR_STATES,
)
from .progress import WorkflowProgressDisplay


class PollingTracker(Protocol):

    def sleep(self) -> None: ...


class PollingTrackerImpl(PollingTracker):

    def __init__(self, polling_backoff: int, timeout=None):
        self.polling_backoff = polling_backoff
        self.timeout = timeout
        self.delta = 0.25
        self.total_wait_time = 0

    def sleep(self):
        if self.timeout is not None and self.total_wait_time > self.timeout:
            message = "Timed out while polling Galaxy."
            raise Exception(message)
        self.total_wait_time += self.delta
        time.sleep(self.delta)
        self.delta += self.polling_backoff


def wait_for_invocation_and_jobs(
    ctx,
    invocation_id: str,
    invocation_api: InvocationApi,
    polling_tracker: PollingTracker,
    workflow_progress_display: WorkflowProgressDisplay,
):
    ctx.vlog("Waiting for invocation [%s]" % invocation_id)

    def summarize(invocation_id: str):
        invocation = invocation_api.get_invocation(invocation_id)
        assert invocation
        invocation_jobs = invocation_api.get_invocation_summary(invocation_id)
        return invocation, invocation_jobs

    last_invocation = None
    last_invocation_jobs = None
    last_subworkflow_invocation = None
    last_subworkflow_invocation_jobs = None
    last_exception = None
    error_message: Optional[str] = None

    done_polling = False
    while not done_polling:
        # loop over the main workflow and one subworkflow each iteration for display,

        # skip the main workflow if it is already tracked as complete - if all steps have been
        # scheduled there are no new subworkflow invocations to track.
        if not workflow_progress_display.workflow_progress.terminal:
            try:
                last_invocation, last_invocation_jobs = summarize(invocation_id)
                workflow_progress_display.handle_invocation(last_invocation, last_invocation_jobs)
            except Exception as e:
                print(e)
                last_exception = e

            error_message = workflow_in_error_message(
                ctx,
                invocation_id,
                last_exception,
                last_invocation,
                last_invocation_jobs,
            )
            if error_message:
                final_invocation_state = "new" if not last_invocation else last_invocation["state"]
                job_state = summary_job_state(last_invocation_jobs)
                return final_invocation_state, job_state, error_message

        assert last_invocation  # if we got here... the first check has passed and we have an invocation

        # grab a subworkflow that isn't complete and check it, also register its subworkflow
        # invocations so we catch all the children and children of children...
        if not workflow_progress_display.all_subworkflows_complete():
            try:
                a_subworkflow_invocation_id = workflow_progress_display.an_incomplete_subworkflow_id()
                last_subworkflow_invocation, last_subworkflow_invocation_jobs = summarize(a_subworkflow_invocation_id)
                workflow_progress_display.handle_subworkflow_invocation(
                    last_subworkflow_invocation, last_subworkflow_invocation_jobs
                )
            except Exception as e:
                last_exception = e

            error_message = workflow_in_error_message(
                ctx,
                invocation_id,
                last_exception,
                last_subworkflow_invocation,
                last_subworkflow_invocation_jobs,
            )
            if error_message:
                final_invocation_state = (
                    "new" if not last_subworkflow_invocation else last_subworkflow_invocation["state"]
                )
                job_state = summary_job_state(last_subworkflow_invocation_jobs)
                return final_invocation_state, job_state, error_message

        done_polling = (
            workflow_progress_display.workflow_progress.terminal
            and workflow_progress_display.all_subworkflows_complete()
        )
        if not done_polling:
            polling_tracker.sleep()

    ctx.vlog(f"The final state of all jobs and subworkflow invocations for invocation [{invocation_id}] is 'ok'")
    job_state = summary_job_state(last_invocation_jobs)
    assert last_invocation
    return last_invocation["state"], job_state, error_message


def workflow_in_error_message(
    ctx, invocation_id, last_exception, last_invocation, last_invocation_jobs
) -> Optional[str]:
    """Return an error message if workflow is in an error state."""

    invocation_state = "new" if not last_invocation else last_invocation["state"]
    job_state = summary_job_state(last_invocation_jobs)

    error_message = None
    if last_exception:
        ctx.vlog(f"Problem waiting on invocation: {str(last_exception)}")
        error_message = f"Final state of invocation {invocation_id} is [{invocation_state}]"

    if invocation_state_terminal(invocation_state) and invocation_state != "scheduled":
        msg = f"Failed to run workflow, invocation ended in [{invocation_state}] state."
        ctx.vlog(msg)
        error_message = msg if not error_message else f"{error_message}. {msg}"

    if job_state in JOB_ERROR_STATES:
        msg = f"Failed to run workflow, at least one job is in [{job_state}] state."
        ctx.vlog(msg)
        error_message = msg if not error_message else f"{error_message}. {msg}"

    return error_message


# we're still mocking out the old history state by just picking out a random
# job state of interest. Seems like we should drop this.
def summary_job_state(job_states_summary: Optional[InvocationJobsSummary]):
    states = (job_states_summary or {"states": {}}).get("states", {}).copy()
    states.pop("ok", None)
    states.pop("skipped", None)
    if states:
        return next(iter(states.keys()))
    else:
        return "ok"


def subworkflow_invocation_ids(invocation_api: InvocationApi, invocation_id: str) -> List[str]:
    invocation = invocation_api.get_invocation(invocation_id)
    subworkflow_invocation_ids = []
    for step in invocation["steps"]:
        subworkflow_invocation_id = step.get("subworkflow_invocation_id")
        if subworkflow_invocation_id:
            subworkflow_invocation_ids.append(subworkflow_invocation_id)
    return subworkflow_invocation_ids
