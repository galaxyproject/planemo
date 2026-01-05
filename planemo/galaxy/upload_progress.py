"""Progress display for Galaxy data uploads."""

from io import StringIO
from typing import (
    Any,
    Dict,
    List,
    Optional,
    Set,
    Union,
)

from rich.console import (
    Console,
    Group,
)
from rich.live import Live
from rich.panel import Panel
from rich.progress import (
    BarColumn,
    Progress,
    TaskID,
    TaskProgressColumn,
    TextColumn,
)

from .invocations.api import JOB_ERROR_STATES
from .invocations.progress_display import DisplayConfiguration

# Job state constants for uploads (same as workflow jobs)
UPLOAD_JOB_TERMINAL_STATES = JOB_ERROR_STATES + ["ok", "skipped", "paused"]


class UploadJobsSummary:
    """Summary of upload job states."""

    def __init__(self, states: Optional[Dict[str, int]] = None):
        self.states: Dict[str, int] = states or {}


def count_states(job_summary: Optional[UploadJobsSummary], query_states: List[str]) -> int:
    """Count jobs in specific states."""
    count = 0
    if job_summary and job_summary.states:
        for state in query_states:
            count += job_summary.states.get(state, 0)
    return count


def job_count(job_summary: Optional[UploadJobsSummary]) -> int:
    """Count total jobs across all states."""
    count = 0
    if job_summary and job_summary.states:
        for state_count in job_summary.states.values():
            if state_count:
                count += state_count
    return count


def ok_count(job_summary: UploadJobsSummary) -> int:
    """Count successfully completed jobs."""
    return count_states(job_summary, ["ok", "skipped"])


def error_count(job_summary: UploadJobsSummary) -> int:
    """Count jobs in error states."""
    return count_states(job_summary, JOB_ERROR_STATES)


class UploadProgress(Progress):
    """Progress tracker for upload jobs."""

    _uploads_task: TaskID

    def __init__(self, display: DisplayConfiguration):
        self.job_count: int = 0
        self.jobs_completed: int = 0
        self.num_ok: int = 0
        self.num_new: int = 0
        self.num_queued: int = 0
        self.num_running: int = 0
        self.num_errors: int = 0
        self.num_paused: int = 0
        self.printed_job_errors: Set[str] = set()  # Track printed job errors by job ID

        self.display = display
        self.uploads_color: str = self.display.style_initializing

        # Create progress bar with same column layout as WorkflowProgress
        bar_column = BarColumn(
            style=self.display.style_bar_back,
            finished_style=self.display.style_bar_finished,
            complete_style=self.display.style_bar_complete,
        )

        super().__init__(
            TextColumn("[progress.description]{task.description}"),
            TextColumn(display.divider),
            bar_column,
            TextColumn(display.divider),
            TaskProgressColumn(f"[{self.display.style_percent}]{{task.percentage:>3.0f}}%"),
            TextColumn(display.divider),
            TextColumn(text_format="{task.fields[status]}"),
        )

        # Add the uploads progress bar
        self._uploads_task = self.add_task(f"[{self.uploads_color}]{self.display.label_progress_jobs}", status="")

    @property
    def jobs_terminal_count(self) -> int:
        """Count of jobs in terminal states (completed, errors, paused)."""
        return self.jobs_completed + self.num_paused

    @property
    def terminal(self) -> bool:
        """Check if all upload jobs are in terminal state."""
        return self.job_count > 0 and self.job_count == self.jobs_terminal_count

    def handle_upload_jobs(self, job_summary: UploadJobsSummary):
        """Update progress based on upload job states.

        Args:
            job_summary: Summary of current job states
        """
        # Update state counters
        self.job_count = job_count(job_summary)
        self.num_new = count_states(job_summary, ["new", "upload"])
        self.num_queued = count_states(job_summary, ["queued", "waiting"])
        self.num_running = count_states(job_summary, ["running"])
        self.num_errors = error_count(job_summary)
        self.num_ok = ok_count(job_summary)
        self.jobs_completed = self.num_ok + self.num_errors
        self.num_paused = count_states(job_summary, ["paused"])

        # Update color based on state
        if self.num_errors > 0:
            self.uploads_color = self.display.style_error
        elif self.job_count > 0:
            self.uploads_color = self.display.style_ok
        else:
            self.uploads_color = self.display.style_initializing

        # Format status message
        jobs_status = ""
        jobs_total: Optional[int] = self.job_count
        jobs_completed: Optional[int] = self.jobs_completed

        if self.job_count == 0:
            jobs_completed = None
            jobs_total = None
        elif self.job_count > 0:
            jobs_status = f"{self.jobs_completed}/{self.job_count} terminal"

        # Update the progress bar
        self.update(
            self._uploads_task,
            total=jobs_total,
            completed=jobs_completed,
            description=f"[{self.uploads_color}]{self.display.label_progress_jobs}",
            status=jobs_status,
        )

    def _job_states_console_line(self) -> str:
        """Format a console line showing breakdown of job states with icons."""
        output = StringIO()
        if self.num_ok > 0:
            output.write(f"{self.display.icon_state_ok} {self.num_ok} {self.display.divider} ")
        if self.num_errors > 0:
            output.write(f"{self.display.icon_state_errors} {self.num_errors} {self.display.divider} ")
        if self.num_new > 0:
            output.write(f"{self.display.icon_state_new} {self.num_new} {self.display.divider} ")
        if self.num_queued > 0:
            output.write(f"{self.display.icon_state_queued} {self.num_queued} {self.display.divider} ")
        if self.num_running > 0:
            output.write(f"{self.display.icon_state_running} {self.num_running} {self.display.divider} ")
        if self.num_paused > 0:
            output.write(f"{self.display.icon_state_paused} {self.num_paused} {self.display.divider} ")

        result = output.getvalue().rstrip(f" {self.display.divider} ")
        output.close()
        return f"[{self.uploads_color}]{self.display.label_job_states_prefix} [reset]{self.display.divider} {result}"

    def print_job_errors_once(self, gi, upload_progress_display: "UploadProgressDisplay"):
        """Print job errors only if they haven't been printed before, tracking by job ID.

        Args:
            gi: GalaxyInstance for querying job details
            upload_progress_display: Display instance for console output
        """
        # Early exit if no errors detected
        if self.num_errors == 0:
            return

        try:
            new_error_lines = []

            # Get all tracked upload jobs and check for errors
            for job_id in upload_progress_display._upload_job_ids:
                if job_id in self.printed_job_errors:
                    continue

                try:
                    # Get job details
                    job_details = gi.jobs.show_job(job_id, full_details=True)
                    job_state = job_details.get("state", "")

                    # If this job is in error state and we haven't printed it yet
                    if job_state in JOB_ERROR_STATES:
                        self.printed_job_errors.add(job_id)
                        error_lines = self._format_job_error_details(job_id, job_details)
                        new_error_lines.extend(error_lines)
                except Exception:
                    # Silently ignore job query failures
                    pass

            # Print any new errors found
            if new_error_lines:
                upload_progress_display.console.print("\n".join(new_error_lines))

        except Exception as e:
            error_msg = f"Failed to collect upload job error details: {e}"
            upload_progress_display.console.print(error_msg)

    def _format_job_error_details(self, job_id: str, job_details: Dict[str, Any]) -> List[str]:
        """Format error details for a single failed upload job.

        Args:
            job_id: Job ID
            job_details: Job details from Galaxy API

        Returns:
            List of formatted error message lines
        """
        error_lines = []
        exit_code = job_details.get("exit_code")
        stderr = job_details.get("stderr", "").strip()
        stdout = job_details.get("stdout", "").strip()
        tool_id = job_details.get("tool_id")

        error_lines.append(f"Failed upload job {job_id}:")

        if tool_id:
            error_lines.append(f"  Tool: {tool_id}")

        if exit_code is not None:
            error_lines.append(f"  Exit code: {exit_code}")

        if stderr:
            error_lines.append(f"  Stderr: {stderr[:500]}")  # Limit stderr output

        if stdout:
            error_lines.append(f"  Stdout: {stdout[:500]}")  # Limit stdout output

        return error_lines


class UploadProgressDisplay(Live):
    """Live display for upload progress with Rich panel."""

    def __init__(
        self,
        history_id: str,
        display_configuration: Optional[DisplayConfiguration] = None,
        galaxy_url: Optional[str] = None,
    ):
        """Initialize upload progress display.

        Args:
            history_id: Galaxy history ID where files are being uploaded
            display_configuration: Optional display configuration for styling
            galaxy_url: Optional Galaxy server URL for display
        """
        self.history_id = history_id
        self.galaxy_url = galaxy_url
        display = display_configuration or DisplayConfiguration()
        self.display = display
        self.upload_progress = UploadProgress(display)
        self.console = Console()
        self._upload_job_ids: Set[str] = set()

        # Initialize with auto_refresh=False for manual control
        super().__init__(self._panel(), console=self.console, auto_refresh=False)

    def get_history_ui_link(self) -> Optional[str]:
        """Get the Galaxy UI link for the history.

        Returns:
            URL to the history in Galaxy UI, or None if galaxy_url not set
        """
        if self.galaxy_url:
            return f"{self.galaxy_url}/histories/view?id={self.history_id}"
        else:
            return None

    def _panel(self) -> Panel:
        """Create Rich panel with progress display.

        Returns:
            Rich Panel containing the upload progress
        """
        # Format title with clickable link like "Uploading to History <hist123abc>"
        history_link = self.get_history_ui_link()
        title = f"[{self.display.style_header}]Uploading to History <[link={history_link}]{self.history_id}[/link]>"

        # Create group with progress bar and optional job state breakdown
        renderables: List[Union[UploadProgress, str]] = [self.upload_progress]

        # Add job state breakdown if configured
        if self.display.include_job_state_breakdown and self.upload_progress.job_count > 0:
            renderables.append(self.upload_progress._job_states_console_line())

        content = Group(*renderables) if len(renderables) > 1 else self.upload_progress

        return Panel(content, title=title, expand=True)

    def update_jobs(self, upload_jobs: List[Dict[str, Any]], job_summary: UploadJobsSummary):
        """Update progress display with current upload job states.

        Args:
            upload_jobs: List of upload job dictionaries
            job_summary: Summary of job states
        """
        # Track job IDs for error reporting
        for job in upload_jobs:
            if "id" in job:
                self._upload_job_ids.add(job["id"])

        # Update progress
        self.upload_progress.handle_upload_jobs(job_summary)

        # Refresh the display
        self.update(self._panel())
        self.refresh()


def _aggregate_job_states(upload_jobs: List[Dict[str, Any]], gi) -> UploadJobsSummary:
    """Poll all upload jobs and aggregate their states.

    Args:
        upload_jobs: List of upload job dictionaries with 'id' field
        gi: GalaxyInstance for querying job states

    Returns:
        UploadJobsSummary with aggregated state counts
    """
    states: Dict[str, int] = {}

    for upload_job in upload_jobs:
        job_id = upload_job.get("id")
        if not job_id:
            continue

        try:
            job = gi.jobs.show_job(job_id)
            state = job.get("state", "unknown")
            states[state] = states.get(state, 0) + 1
        except Exception:
            # If we can't get job state, count it as unknown
            states["unknown"] = states.get("unknown", 0) + 1

    return UploadJobsSummary(states)
