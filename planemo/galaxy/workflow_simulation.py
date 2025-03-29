"""Simulate Galaxy workflows running on a server for testing purposes."""

import yaml
from collections import deque
from typing import (
    List,
    Optional,
)
from uuid import uuid4

from .workflow_progress import (
    Invocation as InvocationResponse,
    InvocationJobsSummary,
    WorkflowProgressDisplay,
)

class Ticks:
    after: int = 0

    @property
    def active(self):
        return self.after <= 0

    def tick(self) -> None:
        if self.active:
            self.tick_when_active()
        else:
            self.after -= 1

    def tick_when_active(self) -> None: ...



class StateWithDuration(Ticks):
    after = 0

    def __init__(self, state: str, duration: int):
        self.state = state
        self.duration = duration

    def tick_when_active(self) -> None:
        self.duration -= 1


class HasState(Ticks):
    final: Optional[str]

    def __init__(self, after: int, states: List[StateWithDuration]):
        self.after = after or 0
        self.states = deque(states)
        self.final_state = None

    def tick_when_active(self) -> None:
        if self.final_state is not None:
            return

        next_state = self.states.popleft()
        next_state.tick()
        if next_state.duration == 0 and not self.states:
            self.final_state = next_state.state
        elif next_state.duration != 0:
            self.states.appendleft(next_state)
        # else: next state will be state during next tick

    @property
    def state(self):
        if self.final_state is not None:
            return self.final_state
        else:
            return self.states[0].state


Job = HasState



class InvocationStep(HasState):
    invocation: Optional["Invocation"]
    jobs: Optional[List[Job]]

    def __init__(self, jobs: List[Job], invocation: Optional["Invocation"], after: int, states: List[StateWithDuration]):
        super().__init__(after, states)
        self.jobs = jobs
        self.invocation = invocation

    def tick_when_active(self) -> None:
        super().tick_when_active()
        if self.jobs:
            for job in self.jobs:
                job.tick()
        if self.invocation:
            self.invocation.tick()
    
    @property
    def active_jobs(self) -> List[Job]:
        return [j for j in (self.jobs or []) if j.active]


class Invocation(HasState):

    def __init__(self, steps: List[InvocationStep], after: int, states: List[StateWithDuration]):
        self.id = str(uuid4())[:10]
        self.steps = steps
        super().__init__(after, states)

    def tick_when_active(self) -> None:
        super().tick_when_active()
        for step in self.steps:
            step.tick()

    @property
    def active_steps(self) -> List[InvocationStep]:
        return [s for s in self.steps if s.active]

    def get_invocation(self) -> InvocationResponse:
        steps = []
        for step in self.active_steps:
            api_step = {
                "state": step.state,
            }
            if step.invocation:
                api_step["subworkflow_invocation_id"] = step.invocation.id

            steps.append(api_step)
        return {
            "id": self.id,
            "state": self.state,
            "steps": steps,
        }

    def get_subworkflow_invocation(self, subworkflow_invocation_id: str) -> "Invocation":
        for step in self.steps:
            if step.invocation and step.invocation.id == subworkflow_invocation_id:
                return step.invocation
        raise Exception(f"Unknown subworkflow invocation id ({subworkflow_invocation_id})")

    def get_subworkflow_invocation_by_step_index(self, index: int) -> "Invocation":
        return self.steps[index].invocation

    def get_jobs_summary(self) -> InvocationJobsSummary:
        jobs = []
        for step in self.active_steps:
            for job in step.active_jobs:
                api_job = {
                    "state": job.state,
                }
                jobs.append(api_job)
        by_state = {}
        for job in jobs:
            state = job["state"]
            if state not in by_state:
                by_state[state] = 0
            by_state[state] += 1
        return {"states": by_state}


def parse_workflow_simulation_from_string(workflow_simulation: str) -> Invocation:
    return parse_workflow_simulation(yaml.safe_load(workflow_simulation))


def parse_workflow_simulation(workflow_simulation: dict) -> Invocation:
    return parse_workflow_simulation_invocation(workflow_simulation)


def parse_workflow_simulation_job(workflow_simulation_job: dict) -> Job:
    states = parse_states_from(workflow_simulation_job)
    after = parse_after_from(workflow_simulation_job)
    return Job(after, states)


def parse_workflow_simulation_invocation_step(workflow_simulation_invocation_step: dict) -> InvocationStep:
    states = parse_states_from(workflow_simulation_invocation_step)
    after = parse_after_from(workflow_simulation_invocation_step)
    if "invocation" in workflow_simulation_invocation_step:
        invocation = parse_workflow_simulation_invocation(workflow_simulation_invocation_step["invocation"])
    else:
        invocation = None
    jobs = None
    if "jobs" in workflow_simulation_invocation_step:
        jobs = []
        for job in workflow_simulation_invocation_step.get("jobs"):
            jobs.append(parse_workflow_simulation_job(job))
    return InvocationStep(jobs, invocation, after, states)


def parse_workflow_simulation_invocation(workflow_simulation_invocation: dict) -> Invocation:
    states = parse_states_from(workflow_simulation_invocation)
    after = parse_after_from(workflow_simulation_invocation)
    steps = []
    for step in workflow_simulation_invocation.get("steps"):
        steps.append(parse_workflow_simulation_invocation_step(step))

    return Invocation(steps, after, states)


def parse_after_from(simulation_object: dict) -> int:
    return simulation_object.get("after", 0)


def parse_states_from(simulation_object: dict) -> List[StateWithDuration]:
    if "states" in simulation_object:
        states = simulation_object["states"]
        states_with_duration = []
        for state in states:
            if ":" in state:
                state, duration_str = state.split(":", 1)
                duration = int(duration_str)
                state_with_duration = StateWithDuration(state, duration)
            else:
                state_with_duration = StateWithDuration(state, 1)
            states_with_duration.append(state_with_duration)
        return states_with_duration
    else:
        state = simulation_object["state"]
        return [StateWithDuration(state, 1)]


def simple_display(invocation: Invocation, display: WorkflowProgressDisplay):
    do_tick = True
    while do_tick:
        last_invocation = invocation.get_invocation()
        last_invocation_jobs_summary = invocation.get_jobs_summary()
        display.handle_invocation(last_invocation, last_invocation_jobs_summary)
        import time

        time.sleep(1)
        invocation.tick()
