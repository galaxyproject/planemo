from time import sleep
from typing import (
    List,
    Optional,
)

from planemo.galaxy.invocations.api import Invocation as InvocationResponse
from planemo.galaxy.invocations.api import (
    InvocationApi,
    InvocationJobsSummary,
    Job,
)
from planemo.galaxy.invocations.polling import (
    PollingTracker,
    wait_for_invocation_and_jobs,
)
from planemo.galaxy.invocations.progress import WorkflowProgressDisplay
from planemo.galaxy.invocations.progress_display import DisplayConfiguration
from planemo.galaxy.invocations.simulations import (
    Invocation,
    parse_workflow_simulation_from_string,
)
from .test_utils import create_test_context
from .test_workflow_simulation import (
    SCENARIO_1,
    SCENARIO_MULTIPLE_OK_SUBWORKFLOWS,
    SCENARIO_NESTED_SUBWORKFLOWS,
    SCENARIO_SUBWORKFLOW_WITH_FAILED_JOBS,
)

SLEEP = 0


class MockPollingTracker(PollingTracker):
    def __init__(self, simulation: Invocation):
        self._simulation = simulation

    def sleep(self) -> None:
        self._simulation.tick()
        if SLEEP > 0:
            sleep(SLEEP)


def test_polling_scenario_1():
    final_invocation_state, job_state, error_message = run_workflow_simulation(SCENARIO_1, fail_fast=True)
    assert final_invocation_state == "ready"  # early job error and fail fast, invocation doesn't advance to scheduled
    assert job_state == "failed"
    assert error_message
    assert "failed" in error_message


def test_polling_scenario_three_ok_subworkflows():
    final_invocation_state, job_state, error_message = run_workflow_simulation(SCENARIO_MULTIPLE_OK_SUBWORKFLOWS)
    assert final_invocation_state == "scheduled"
    assert job_state == "ok"
    assert not error_message


def test_polling_scenario_nested_subworkflows():
    final_invocation_state, job_state, error_message = run_workflow_simulation(SCENARIO_NESTED_SUBWORKFLOWS)
    assert final_invocation_state == "scheduled"
    assert job_state == "ok"
    assert not error_message


def test_polling_without_display():
    simulation = parse_workflow_simulation_from_string(SCENARIO_1)
    invocation_id = simulation.id
    invocation_api = SimulatedApi(simulation)
    polling_tracker = MockPollingTracker(simulation)
    ctx = create_test_context()
    # using this without setting up the display context seems to use all the tracking
    # without printing to the console. Hides a lot of bad design mixing presentation and
    # tracking logic being too mixed together.
    display = WorkflowProgressDisplay(invocation_id)
    final_invocation_state, job_state, error_message = wait_for_invocation_and_jobs(
        ctx,
        invocation_id,
        invocation_api,
        polling_tracker,
        display,
        fail_fast=True,
    )
    assert final_invocation_state == "ready"
    assert job_state == "failed"
    assert error_message
    assert "failed" in error_message


def test_polling_with_compact_display():
    display_configuration = DisplayConfiguration(
        include_nested_subworkflows=False,
        include_job_state_breakdown=False,
    )
    final_invocation_state, job_state, error_message = run_workflow_simulation(
        SCENARIO_NESTED_SUBWORKFLOWS, display_configuration
    )
    assert final_invocation_state == "scheduled"
    assert job_state == "ok"
    assert not error_message


def test_polling_without_invocation_as_full_subpanel():
    display_configuration = DisplayConfiguration(
        include_nested_subworkflows=True,
        include_job_state_breakdown=True,
        subworkflows_as_panel=False,
    )
    final_invocation_state, job_state, error_message = run_workflow_simulation(
        SCENARIO_NESTED_SUBWORKFLOWS, display_configuration
    )
    assert final_invocation_state == "scheduled"
    assert job_state == "ok"
    assert not error_message


def test_fail_fast_enabled_with_job_failure():
    """Test that fail_fast=True returns error when a job fails."""
    final_invocation_state, job_state, error_message = run_workflow_simulation(SCENARIO_1, fail_fast=True)
    # Invocation should still be scheduled (workflow scheduling succeeded)
    assert final_invocation_state == "ready"
    assert job_state == "failed"
    # fail_fast should detect the failed job and return error message
    assert error_message
    assert "Failed to run workflow, at least one job is in [failed] state." in error_message


def test_fail_fast_disabled_with_job_failure():
    """Test that fail_fast=False does not report job failures as errors."""
    final_invocation_state, job_state, error_message = run_workflow_simulation(SCENARIO_1, fail_fast=False)
    # Invocation should be scheduled (workflow scheduling succeeded)
    assert final_invocation_state == "scheduled"
    assert job_state == "failed"
    # Without fail_fast, job failures shouldn't cause error messages
    # (unless invocation itself fails, which it doesn't in this case)
    assert error_message is None


def test_fail_fast_enabled_with_successful_workflow():
    """Test that fail_fast=True works normally when no jobs fail."""
    final_invocation_state, job_state, error_message = run_workflow_simulation(
        SCENARIO_MULTIPLE_OK_SUBWORKFLOWS, fail_fast=True
    )
    assert final_invocation_state == "scheduled"
    assert job_state == "ok"
    assert not error_message


def test_fail_fast_enabled_with_subworkflow_job_failure():
    """Test that fail_fast=True terminates when encountering jobs that are errored inside a subworkflow invocation."""
    final_invocation_state, job_state, error_message = run_workflow_simulation(
        SCENARIO_SUBWORKFLOW_WITH_FAILED_JOBS, fail_fast=True
    )
    # Invocation is ready to schedule more steps, yet the polling should terminate
    assert final_invocation_state == "ready"
    assert job_state == "error"
    # fail_fast should detect the failed job in the subworkflow and return error message
    assert error_message
    assert "Failed to run workflow, at least one job is in [error] state." in error_message


def run_workflow_simulation(
    yaml_str: str, display_configuration: Optional[DisplayConfiguration] = None, fail_fast: bool = False
):
    simulation = parse_workflow_simulation_from_string(yaml_str)
    invocation_id = simulation.id
    invocation_api = SimulatedApi(simulation)
    polling_tracker = MockPollingTracker(simulation)
    ctx = create_test_context()
    with WorkflowProgressDisplay(invocation_id, display_configuration=display_configuration) as display:
        return wait_for_invocation_and_jobs(
            ctx,
            invocation_id,
            invocation_api,
            polling_tracker,
            display,
            fail_fast=fail_fast,
        )


class MockGalaxyInstance:
    """Mock Galaxy instance for testing job detail retrieval."""

    def __init__(self):
        self.jobs = MockJobsApi()
        self.invocations = MockInvocationsApi()


class MockJobsApi:
    """Mock jobs API for testing."""

    def show_job(self, job_id, full_details=False):
        """Return mock job details with exit code and stderr."""
        return {
            "id": job_id,
            "state": "failed",
            "exit_code": 1,
            "stderr": f"Error: Mock job {job_id} failed with exit code 1\nAdditional error details here",
            "stdout": f"Mock job {job_id} output",
            "command_line": f"python /mock/tool/script.py --input data.txt --output {job_id}_output.txt",
            "tool_id": f"mock_tool_{job_id.split('_')[-1]}",
        }

    def get_jobs(self, invocation_id=None, state=None):
        """Return mock list of jobs with specified state for invocation."""
        if state == "error" and invocation_id:
            return [
                {"id": "failed_job_1", "state": "error"},
                {"id": "failed_job_2", "state": "error"},
                {"id": "failed_job_3", "state": "error"},
                {"id": "failed_job_4", "state": "error"},
            ]
        return []


class MockInvocationsApi:
    """Mock invocations API for testing."""

    def show_invocation_step(self, invocation_id, step_id):
        """Return mock invocation step details."""
        return {"id": step_id, "jobs": [{"id": f"job_{step_id}", "state": "failed"}]}


class SimulatedApi(InvocationApi):
    _simulation: Invocation

    def __init__(self, invocation: Invocation):
        self._simulation = invocation
        self._user_gi = MockGalaxyInstance()

    def get_invocation(self, invocation_id: str) -> InvocationResponse:
        invocation = self._simulation.get_invocation_by_id(invocation_id)
        assert invocation, f"Simulation has no invocation_id {invocation_id}"
        return invocation.get_api_invocation()

    def get_invocation_summary(self, invocation_id: str, state: Optional[str] = None) -> InvocationJobsSummary:
        invocation = self._simulation.get_invocation_by_id(invocation_id)
        assert invocation, f"Simulation has no invocation_id {invocation_id}"
        return invocation.get_api_jobs_summary()

    def get_invocation_jobs(self, invocation_id: str, state: Optional[str] = None) -> List[Job]:
        """Return mock list of jobs for invocation, filtered by state if provided."""
        if state == "error":
            return [
                {"id": "failed_job_1", "state": "error"},
                {"id": "failed_job_2", "state": "error"},
                {"id": "failed_job_3", "state": "error"},
                {"id": "failed_job_4", "state": "error"},
            ]
        return []

    def get_job(self, job_id: str, full_details: bool = False) -> Job:
        """Return mock job details."""
        return {
            "id": job_id,
            "state": "failed",
            "exit_code": 1,
            "stderr": f"Error: Mock job {job_id} failed with exit code 1\nAdditional error details here",
            "stdout": f"Mock job {job_id} output",
            "command_line": f"python /mock/tool/script.py --input data.txt --output {job_id}_output.txt",
            "tool_id": f"mock_tool_{job_id.split('_')[-1]}",
        }
