import random
from typing import (
    Dict,
    List,
    Optional,
    Set,
)

from io import StringIO
from rich.console import Group
from rich.live import Live
from rich.panel import Panel
from rich.console import Console
from rich.progress import (
    BarColumn,
    Progress,
    TextColumn,
    TaskProgressColumn,
)
from rich.theme import Theme
from typing_extensions import TypedDict

# uses a bit more space but has better visual separation between subworkflows and workflows
DISPLAY_SUBWORKFLOW_AS_PANEL = True
DISPLAY_INCLUDE_JOB_STATE_BREAKDOWN = True

# Types for various invocation responses
class InvocationStep(TypedDict, total=False):
    state: Optional[str]
    subworkflow_invocation_id: Optional[str]


class Invocation(TypedDict, total=False):
    id: str
    state: str
    steps: List[InvocationStep]


class InvocationJobsSummary(TypedDict, total=False):
    states: Dict[str, int]


class WorkflowProgress(Progress):
    invocation_state: str = "new"
    step_count: Optional[int] = None
    job_count: int = 0
    steps_color: str = "cyan"
    jobs_color: str = "cyan"
    step_states: Dict = {}
    num_ok: int = 0
    num_new: int = 0
    num_queued: int = 0
    num_running: int = 0
    num_errors: int = 0
    num_paused: int = 0

    def __init__(self):
        custom_theme = Theme(
            {
                "bar.complete": "green",
            }
        )
        console = Console(theme=custom_theme)
        super().__init__(
            TextColumn("[progress.description]{task.description}"),
            TextColumn("◆"),
            BarColumn(),
            TextColumn("◆"),
            TaskProgressColumn(),
            TextColumn("◆"),
            TextColumn(text_format="{task.fields[status]}"),
            console=console,
        )
        self.add_bars()

    @property
    def invocation_scheduling_terminal(self):
        return invocation_state_terminal(self.invocation_state)

    @property
    def jobs_terminal(self):
        return self.job_count is not None and self.job_count == self.jobs_terminal_count

    def handle_invocation(self, invocation: Invocation, job_state_summary: InvocationJobsSummary):
        self.invocation_state = invocation.get("state") or "new"
        self.step_count = len(invocation.get("steps") or []) or None
        self.step_states = step_states(invocation)

        steps_completed = None

        steps_status = ""
        if self.step_count is None:
            steps_status = "Loading steps."
            self.steps_color = "cyan"
        elif self.invocation_state == "cancelled":
            steps_status = "Invocation cancelled"
            self.steps_color = "red"
        elif self.invocation_state == "failed":
            steps_status = "Invocation failed"
            self.steps_color = "red"
        else:
            num_scheduled = self.step_states.get("scheduled") or 0
            if num_scheduled > 0:
                self.steps_color = "green"
            else:
                self.steps_color = "cyan"
            steps_completed = num_scheduled
            steps_status = f"{num_scheduled}/{self.step_count} scheduled"

        jobs_status = ""
        self.job_count = job_count(job_state_summary)
        self.num_new = count_states(job_state_summary, ["new"])
        self.num_queued = count_states(job_state_summary, ["queued", "waiting"])
        self.num_running = count_states(job_state_summary, ["running"])
        self.num_errors = error_count(job_state_summary)
        self.num_ok = ok_count(job_state_summary)
        self.jobs_completed = self.num_ok + self.num_errors
        self.num_paused = count_states(job_state_summary, ["paused"])
        self.jobs_terminal_count = self.jobs_completed + self.num_paused
        jobs_total = self.job_count
        if self.num_errors > 0:
            self.jobs_color = "red"
        elif self.job_count > 0:
            self.jobs_color = "green"
        else:
            self.jobs_color = "cyan"
            self.jobs_completed = None
            jobs_total = None
        if self.job_count > 0:
            jobs_status = f"{self.jobs_completed}/{self.job_count} terminal"
        self.update(
            self._steps_task,
            total=self.step_count,
            completed=steps_completed,
            description=f"[{self.steps_color}]Steps",
            status=steps_status,
        )
        self.update(
            self._jobs_task,
            total=jobs_total,
            completed=self.jobs_completed,
            description=f"[{self.jobs_color}]Jobs",
            status=jobs_status,
        )

    def _job_states_console_line(self):
        output = StringIO()
        if self.num_ok > 0:
            output.write(f"🟢 {self.num_ok} ◆ ")
        if self.num_errors > 0:
            output.write(f"🔴 {self.num_errors} ◆ ")
        if self.num_new > 0:
            output.write(f"🆕 {self.num_new} ◆ ")
        if self.num_queued > 0:
            output.write(f"⏳ {self.num_queued} ◆ ")
        if self.num_running > 0:
            output.write(f"👟 {self.num_running} ◆ ")
        if self.num_paused > 0:
            output.write(f"⏸️ {self.num_paused} ◆ ")

        result = output.getvalue().rstrip(" ◆ ")
        output.close()
        return f"[{self.jobs_color}]Job States   [black]◆ {result}"

    def add_bars(self):
        self._steps_task = self.add_task("[green]Steps", status="")
        self._jobs_task = self.add_task("[green]Jobs", status="")


# converted from Galaxy TypeScript (see util.ts next to WorkflowInvocationState.vue)
def count_states(job_summary: Optional[InvocationJobsSummary], query_states: list[str]) -> int:
    count = 0
    states = job_summary.get("states") if job_summary else None
    if states:
        for state in query_states:
            count += states.get(state, 0)
    return count


def job_count(job_summary: Optional[InvocationJobsSummary]) -> int:
    states = job_summary.get("states") if job_summary else None
    count = 0
    if states:
        for state_count in states.values():
            if state_count:
                count += state_count
    return count


def step_states(invocation: Invocation):
    step_states = {}
    steps = invocation.get("steps") or []
    for step in steps:
        if not step:
            continue
        step_state = step.get("state") or "unknown"
        if step_state not in step_states:
            step_states[step_state] = 0
        step_states[step_state] += 1

    return step_states


def ok_count(job_summary: InvocationJobsSummary) -> int:
    return count_states(job_summary, ["ok", "skipped"])


# ERROR_STATES = ["error", "deleted", "failed"]
ERROR_STATES = ["error", "deleted", "failed", "stopped", "stop", "deleting"]


def error_count(job_summary: InvocationJobsSummary) -> int:
    return count_states(job_summary, ERROR_STATES)


def running_count(job_summary: InvocationJobsSummary) -> int:
    return count_states(job_summary, ["running"])


def invocation_state_terminal(state: str):
    return state in ["scheduled", "cancelled", "failed"]


class WorkflowProgressDisplay(Live):
    workflow_progress: WorkflowProgress
    subworkflow_progress: Optional[WorkflowProgress] = None
    subworkflow_invocation_ids_seen: Set[str] = set()
    subworkflow_invocation_ids_completed: Set[str] = set()
    invocation_id: str
    subworkflow_invocation_id: Optional[str] = None

    def __init__(
        self,
        invocation_id: str,
        display_subworkflow_as_panel: bool = DISPLAY_SUBWORKFLOW_AS_PANEL,
        display_include_job_state_breakdown: bool = DISPLAY_INCLUDE_JOB_STATE_BREAKDOWN,
    ):
        self.invocation_id = invocation_id
        self.workflow_progress = WorkflowProgress()
        self.subworkflow_progress = WorkflowProgress()
        self.display_subworkflow_as_panel = display_subworkflow_as_panel
        self.display_include_job_state_breakdown = display_include_job_state_breakdown
        super().__init__(self._panel())

    def register_subworkflow_invocation_ids(self, ids: List[str]):
        for invocation_id in ids:
            self.subworkflow_invocation_ids_seen.add(invocation_id)

    def complete_subworkflow(self, id: str):
        self.subworkflow_invocation_ids_completed.add(id)

    def an_incomplete_subworkflow_id(self):
        return random.choice(tuple(self.subworkflow_invocation_ids_seen - self.subworkflow_invocation_ids_completed))

    def all_subworkflows_complete(self):
        return len(self.subworkflow_invocation_ids_seen) == len(self.subworkflow_invocation_ids_completed)

    def _panel(self):
        def job_states(workflow_progress):
            if self.display_include_job_state_breakdown:
                return workflow_progress._job_states_console_line()
            else:
                return None

        title = f"[bold]Invocation <{self.invocation_id}>"
        subworkflow_title = None
        if self.subworkflow_invocation_id:
            subworkflow_title = f"[bold]Subworkflow Invocation <{self.subworkflow_invocation_id}>"

        if not self.subworkflow_invocation_id:
            renderable = Group(
                self.workflow_progress,
                job_states(self.workflow_progress),
            )
        elif not self.display_subworkflow_as_panel:
            renderable = as_group(
                self.workflow_progress,
                job_states(self.workflow_progress),
                subworkflow_title,
                self.subworkflow_progress,
                job_states(self.subworkflow_progress),
            )
        else:
            renderable = as_group(
                self.workflow_progress,
                job_states(self.workflow_progress),
                Panel(
                    as_group(self.subworkflow_progress, job_states(self.subworkflow_progress)),
                    title=subworkflow_title,
                ),
            )
        return Panel(renderable, title=title, expand=True)

    def _update_panel(self):
        self.update(self._panel())

    def handle_invocation(self, invocation: Invocation, job_state_summary: InvocationJobsSummary):
        self.workflow_progress.handle_invocation(invocation, job_state_summary)
        self._update_panel()

    def handle_subworkflow_invocation(self, invocation: Invocation, job_state_summary: InvocationJobsSummary):
        self.subworkflow_invocation_id = invocation["id"]
        self.subworkflow_progress.handle_invocation(invocation, job_state_summary)
        self._update_panel()


def as_group(*renderables):
    return Group(*(r for r in renderables if r))
