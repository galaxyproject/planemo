from planemo.galaxy.invocations.simulations import parse_workflow_simulation_from_string

# These scenario simulates a workflow execution timeline with the following structure:
#
# Field explanations:
# - states: Array defining the main workflow invocation's state progression over time
#   Format: [initial_state, intermediate_state:duration, final_state]
#   - "new": Workflow just created
#   - "ready:4": Workflow in ready state for 4 time ticks
#   - "scheduled": Workflow is running/scheduled
#
# - steps: Array of workflow steps that execute in sequence
#   - Each step has a "state" and may contain "jobs" or nested "invocation" (subworkflow)
#   - "after: N": This step starts after N time ticks
#
# - jobs: Array of jobs within a step
#   - states: Job state progression over time
#   - Common states: new -> queued -> running -> ok/failed
#   - Format with duration: "queued:2" means job stays queued for 2 ticks

SCENARIO_1 = """
states: [new, ready:4, scheduled]
steps:
- state: scheduled
  jobs:
  - states: [new, queued:2, running:2, ok]
- after: 2
  state: scheduled
  jobs:
  - states: [new, queued, error]
  - states: [new, queued, ok]
- after: 3
  state: scheduled
  invocation:
    states: [new, ready, scheduled]
    steps:
    - state: scheduled
      jobs:
        - states: [new, queued, ok]
        - states: [new, queued:2, ok]
"""

SCENARIO_MULTIPLE_OK_SUBWORKFLOWS = """
states: [new, ready:4, scheduled]
steps:
- state: scheduled
  jobs:
  - states: [new, queued:2, running:2, ok]
- after: 2
  state: scheduled
  jobs:
  - states: [new, queued, running:2, ok]
  - states: [new, queued, running:4, ok]
- after: 3
  state: scheduled
  invocation:
    states: [new, ready, scheduled]
    steps:
    - state: scheduled
      jobs:
        - states: [new, queued, ok]
        - states: [new, queued:2, running, ok]
- after: 4
  state: scheduled
  invocation:
    states: [new, ready, scheduled]
    steps:
    - state: scheduled
      jobs:
        - states: [new, queued, ok]
        - states: [new, queued:2, ok]
    - after: 3
      state: scheduled
      jobs:
        - states: [new, queued, running, ok]
        - states: [new, queued:3, running:2, ok]
    - after: 4
      state: scheduled
      jobs:
        - states: [new, queued, ok]
        - states: [new, queued:2, ok]
- after: 5
  state: scheduled
  invocation:
    states: [new, ready:2, scheduled]
    steps:
    - state: scheduled
      jobs:
        - states: [new, queued:3, ok]
        - states: [new, queued:1, ok]
    - after: 3
      state: scheduled
      jobs:
        - states: [new, queued:2, running:4, ok]
        - states: [new, queued:3, running:2, ok]
    - after: 4
      state: scheduled
      jobs:
        - states: [new, queued, ok]
        - states: [new, running:5, ok]
"""

SCENARIO_NESTED_SUBWORKFLOWS = """
states: [new, ready:4, scheduled]
steps:
- state: scheduled
  jobs:
  - states: [new, queued:2, running:2, ok]
- after: 2
  state: scheduled
  jobs:
  - states: [new, queued, running:2, ok]
  - states: [new, queued, running:4, ok]
- after: 3
  state: scheduled
  invocation:
    states: [new, ready, scheduled]
    steps:
    - state: scheduled
      jobs:
        - states: [new, queued, ok]
        - states: [new, queued:2, running, ok]
    - after: 3
      states: [new, ready, scheduled]
      invocation:
        states: [new, ready:4, scheduled]
        steps:
        - state: scheduled
          jobs:
          - states: [new, queued, ok]
          - states: [new, queued:2, running, ok]
        - state: scheduled
          jobs:
          - states: [new, queued, ok]
          - states: [new, queued:2, running, ok]
          - states: [new, queued:3, running, ok]
- after: 4
  state: scheduled
  invocation:
    states: [new, ready, scheduled]
    steps:
    - state: scheduled
      jobs:
        - states: [new, queued, ok]
        - states: [new, queued:2, ok]
    - after: 3
      state: scheduled
      jobs:
        - states: [new, queued, running, ok]
        - states: [new, queued:3, running:2, ok]
    - after: 4
      state: scheduled
      jobs:
        - states: [new, queued, ok]
        - states: [new, queued:2, ok]
- after: 5
  state: scheduled
  invocation:
    states: [new, ready:2, scheduled]
    steps:
    - state: scheduled
      jobs:
        - states: [new, queued:3, ok]
        - states: [new, queued:1, ok]
    - after: 3
      state: scheduled
      jobs:
        - states: [new, queued:2, running:4, ok]
        - states: [new, queued:3, running:2, ok]
    - after: 4
      state: scheduled
      jobs:
        - states: [new, queued, ok]
        - states: [new, running:5, ok]
"""


SCENARIO_SUBWORKFLOW_WITH_FAILED_JOBS = """
states: [new, ready:4]
steps:
- state: scheduled
  jobs:
  - states: [new, queued:2, running:2, ok]
- after: 2
  state: scheduled
  jobs:
  - states: [new, queued, running:2, ok]
  - states: [new, queued, running:4, ok]
- after: 3
  state: scheduled
  invocation:
    states: [new, ready]
    steps:
    - state: scheduled
      jobs:
        - states: [new, queued, error]
        - states: [new:2, paused]
"""


def test_parse_scenario_1_invocation_state_evolution():
    invocation = parse_workflow_simulation_from_string(SCENARIO_1)
    invocation_dict = invocation.get_api_invocation()
    assert invocation_dict["state"] == "new"
    invocation.tick()
    invocation_dict = invocation.get_api_invocation()
    assert invocation_dict["state"] == "ready"
    invocation.tick()
    invocation.tick()
    invocation.tick()
    invocation.tick()
    invocation_dict = invocation.get_api_invocation()
    assert invocation_dict["state"] == "scheduled"


def test_parse_scenario_1_invocation_step_states():
    invocation = parse_workflow_simulation_from_string(SCENARIO_1)
    invocation_dict = invocation.get_api_invocation()
    steps = invocation_dict["steps"]
    assert len(steps) == 1

    invocation.tick()
    invocation.tick()

    invocation_dict = invocation.get_api_invocation()
    steps = invocation_dict["steps"]
    assert len(steps) == 2
    assert steps[0]["state"] == "scheduled"
    assert steps[1]["state"] == "scheduled"

    invocation.tick()

    invocation_dict = invocation.get_api_invocation()
    steps = invocation_dict["steps"]
    assert len(steps) == 3
    assert steps[2]["state"] == "scheduled"


def test_parse_scenario_1_invocation_job_states():
    invocation = parse_workflow_simulation_from_string(SCENARIO_1)
    states = invocation.get_api_jobs_summary()["states"]
    assert len(states) == 1
    assert states["new"] == 1

    invocation.tick()

    states = invocation.get_api_jobs_summary()["states"]
    assert len(states) == 1
    assert states["queued"] == 1

    invocation.tick()

    states = invocation.get_api_jobs_summary()["states"]
    assert len(states) == 2
    assert states["queued"] == 1
    assert states["new"] == 2

    invocation.tick()

    states = invocation.get_api_jobs_summary()["states"]
    assert len(states) == 2
    assert states["queued"] == 2
    assert states["running"] == 1

    invocation.tick()

    states = invocation.get_api_jobs_summary()["states"]
    assert len(states) == 3
    assert states["ok"] == 1
    assert states["running"] == 1
    assert states["error"] == 1


def test_parse_scenario_1_subworkflow_invocation_state():
    invocation = parse_workflow_simulation_from_string(SCENARIO_1)

    invocation.tick()
    invocation.tick()
    invocation.tick()

    subworkflow_invocation = invocation.get_subworkflow_invocation_by_step_index(2)
    assert subworkflow_invocation.get_api_invocation()["state"] == "new"

    invocation.tick()

    assert subworkflow_invocation.get_api_invocation()["state"] == "ready"

    invocation.tick()

    assert subworkflow_invocation.get_api_invocation()["state"] == "scheduled"

    states = subworkflow_invocation.get_api_jobs_summary()["states"]
    assert len(states) == 2
    assert states["ok"] == 1
    assert states["queued"] == 1

    invocation.tick()

    states = subworkflow_invocation.get_api_jobs_summary()["states"]
    assert len(states) == 1
    assert states["ok"] == 2
