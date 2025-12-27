"""Tests for upload progress display."""

import time

from planemo.galaxy.invocations.progress_display import DisplayConfiguration
from planemo.galaxy.upload_progress import (
    UploadJobsSummary,
    UploadProgress,
    UploadProgressDisplay,
)

# Test timing (sleep duration for visual tests)
SLEEP = 0.8


def test_upload_progress_typical():
    """Test typical upload progress sequence with state transitions."""
    with UploadProgress(DisplayConfiguration()) as upload_progress:
        # Initial state: 3 new jobs
        upload_progress.handle_upload_jobs(UploadJobsSummary({"new": 3}))
        assert not upload_progress.terminal
        assert upload_progress.job_count == 3
        assert upload_progress.num_new == 3
        assert upload_progress.jobs_completed == 0

        # Some jobs move to queued
        upload_progress.handle_upload_jobs(UploadJobsSummary({"queued": 2, "new": 1}))
        assert not upload_progress.terminal
        assert upload_progress.num_queued == 2
        assert upload_progress.num_new == 1

        # Some jobs start running
        upload_progress.handle_upload_jobs(UploadJobsSummary({"running": 2, "new": 1}))
        assert not upload_progress.terminal
        assert upload_progress.num_running == 2

        # First job completes
        upload_progress.handle_upload_jobs(UploadJobsSummary({"running": 1, "ok": 2}))
        assert not upload_progress.terminal
        assert upload_progress.num_ok == 2
        assert upload_progress.jobs_completed == 2

        # All jobs complete
        upload_progress.handle_upload_jobs(UploadJobsSummary({"ok": 3}))
        assert upload_progress.terminal
        assert upload_progress.num_ok == 3
        assert upload_progress.jobs_completed == 3


def test_upload_progress_terminal_detection():
    """Test various terminal state conditions."""
    progress = UploadProgress(DisplayConfiguration())

    # All ok - terminal
    progress.handle_upload_jobs(UploadJobsSummary({"ok": 5}))
    assert progress.terminal
    assert progress.job_count == 5
    assert progress.jobs_completed == 5

    # Mix of ok and error - terminal
    progress.handle_upload_jobs(UploadJobsSummary({"ok": 3, "error": 2}))
    assert progress.terminal
    assert progress.job_count == 5
    assert progress.jobs_completed == 5

    # Some still running - not terminal
    progress.handle_upload_jobs(UploadJobsSummary({"ok": 3, "running": 2}))
    assert not progress.terminal
    assert progress.job_count == 5
    assert progress.jobs_completed == 3

    # Some still new - not terminal
    progress.handle_upload_jobs(UploadJobsSummary({"ok": 3, "new": 2}))
    assert not progress.terminal

    # Some queued - not terminal
    progress.handle_upload_jobs(UploadJobsSummary({"ok": 3, "queued": 2}))
    assert not progress.terminal

    # All paused (terminal)
    progress.handle_upload_jobs(UploadJobsSummary({"paused": 5}))
    assert progress.terminal

    # Mix of completed and paused (terminal)
    progress.handle_upload_jobs(UploadJobsSummary({"ok": 3, "error": 1, "paused": 1}))
    assert progress.terminal


def test_upload_progress_with_errors():
    """Test upload progress with job failures."""
    progress = UploadProgress(DisplayConfiguration())

    # Some jobs running, one error
    progress.handle_upload_jobs(UploadJobsSummary({"running": 2, "error": 1}))
    assert progress.num_errors == 1
    assert not progress.terminal
    assert progress.uploads_color == progress.display.style_error  # Should be red

    # More jobs complete, error remains
    progress.handle_upload_jobs(UploadJobsSummary({"ok": 2, "error": 1}))
    assert progress.terminal
    assert progress.num_errors == 1
    assert progress.num_ok == 2


def test_upload_progress_state_counters():
    """Test that state counters are updated correctly."""
    progress = UploadProgress(DisplayConfiguration())

    # Test various state combinations
    progress.handle_upload_jobs(
        UploadJobsSummary(
            {
                "new": 2,
                "upload": 1,  # upload is counted as "new"
                "queued": 3,
                "waiting": 1,  # waiting is counted as "queued"
                "running": 2,
                "ok": 4,
                "skipped": 1,  # skipped is counted as "ok"
                "error": 1,
            }
        )
    )

    assert progress.num_new == 3  # new + upload
    assert progress.num_queued == 4  # queued + waiting
    assert progress.num_running == 2
    assert progress.num_ok == 5  # ok + skipped
    assert progress.num_errors == 1
    assert progress.job_count == 15
    assert progress.jobs_completed == 6  # ok + errors


def test_upload_progress_color_coding():
    """Test that progress bar color changes based on state."""
    progress = UploadProgress(DisplayConfiguration())

    # No jobs - initializing color
    progress.handle_upload_jobs(UploadJobsSummary({}))
    assert progress.uploads_color == progress.display.style_initializing

    # Jobs running - ok color
    progress.handle_upload_jobs(UploadJobsSummary({"running": 2}))
    assert progress.uploads_color == progress.display.style_ok

    # Jobs with errors - error color
    progress.handle_upload_jobs(UploadJobsSummary({"running": 1, "error": 1}))
    assert progress.uploads_color == progress.display.style_error

    # All ok - ok color
    progress.handle_upload_jobs(UploadJobsSummary({"ok": 2}))
    assert progress.uploads_color == progress.display.style_ok


def test_upload_progress_display_initialization():
    """Test UploadProgressDisplay initialization."""
    display = UploadProgressDisplay("hist123")
    assert display.history_id == "hist123"
    assert display.galaxy_url is None
    assert display.upload_progress is not None
    assert len(display._upload_job_ids) == 0

    # Test with Galaxy URL
    display = UploadProgressDisplay("hist456", galaxy_url="http://localhost:8080")
    assert display.history_id == "hist456"
    assert display.galaxy_url == "http://localhost:8080"


def test_upload_progress_display_update():
    """Test updating upload progress display."""
    display = UploadProgressDisplay("hist789")

    # Update with upload jobs
    upload_jobs = [
        {"id": "job1"},
        {"id": "job2"},
        {"id": "job3"},
    ]
    job_summary = UploadJobsSummary({"new": 3})

    display.update_jobs(upload_jobs, job_summary)

    # Verify job IDs are tracked
    assert "job1" in display._upload_job_ids
    assert "job2" in display._upload_job_ids
    assert "job3" in display._upload_job_ids

    # Verify progress is updated
    assert display.upload_progress.job_count == 3
    assert display.upload_progress.num_new == 3


def test_upload_progress_display_typical():
    """Visual test of upload progress display with state transitions."""
    with UploadProgressDisplay("hist123abcde") as live:
        upload_jobs = [
            {"id": "upload1"},
            {"id": "upload2"},
            {"id": "upload3"},
        ]

        state_sequences = [
            UploadJobsSummary({"new": 3}),
            UploadJobsSummary({"new": 2, "queued": 1}),
            UploadJobsSummary({"new": 1, "queued": 1, "running": 1}),
            UploadJobsSummary({"queued": 1, "running": 2}),
            UploadJobsSummary({"running": 2, "ok": 1}),
            UploadJobsSummary({"running": 1, "ok": 2}),
            UploadJobsSummary({"ok": 3}),
        ]

        for job_summary in state_sequences:
            live.update_jobs(upload_jobs, job_summary)
            time.sleep(SLEEP)


def test_upload_progress_empty_jobs():
    """Test handling of empty job list."""
    progress = UploadProgress(DisplayConfiguration())

    # Handle empty job summary
    progress.handle_upload_jobs(UploadJobsSummary({}))
    assert progress.job_count == 0
    assert not progress.terminal  # Not terminal when no jobs
    assert progress.jobs_completed == 0


def test_upload_progress_job_error_tracking():
    """Test that job errors are tracked to avoid duplicates."""
    progress = UploadProgress(DisplayConfiguration())

    # Initially no errors printed
    assert len(progress.printed_job_errors) == 0

    # Add error to tracking
    progress.printed_job_errors.add("job123")
    assert "job123" in progress.printed_job_errors

    # Verify we can check if error was already printed
    assert "job123" in progress.printed_job_errors
    assert "job456" not in progress.printed_job_errors
