from time import sleep
from typing import Optional

from planemo.galaxy.invocations.api import Invocation as InvocationResponse
from planemo.galaxy.invocations.api import (
    InvocationApi,
    InvocationJobsSummary,
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
    final_invocation_state, job_state, error_message = run_workflow_simulation(SCENARIO_1)
    assert final_invocation_state == "scheduled"
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
    )
    assert final_invocation_state == "scheduled"
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


def run_workflow_simulation(yaml_str: str, display_configuration: Optional[DisplayConfiguration] = None):
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
        )


class SimulatedApi(InvocationApi):
    _simulation: Invocation

    def __init__(self, invocation: Invocation):
        self._simulation = invocation

    def get_invocation(self, invocation_id: str) -> InvocationResponse:
        invocation = self._simulation.get_invocation_by_id(invocation_id)
        assert invocation, f"Simulation has no invocation_id {invocation_id}"
        return invocation.get_api_invocation()

    def get_invocation_summary(self, invocation_id: str) -> InvocationJobsSummary:
        invocation = self._simulation.get_invocation_by_id(invocation_id)
        assert invocation, f"Simulation has no invocation_id {invocation_id}"
        return invocation.get_api_jobs_summary()
