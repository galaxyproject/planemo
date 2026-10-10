"""Tests for retaining an engine after a Planemo test run."""

import contextlib
from types import SimpleNamespace
from unittest.mock import patch

import click
import pytest

from planemo.engine.interface import BaseEngine
from planemo.engine.test import test_runnables as run_tests
from .test_utils import create_test_context


class _UnsupportedEngine(BaseEngine):
    def _run(self, runnables, job_path, output_collectors=None, test_timeout=None):
        raise AssertionError("an unsupported engine must reject --serve before running tests")


class _RecordingEngine:
    can_serve_test_results = True

    def __init__(self, events, test_data):
        self.events = events
        self.test_data = test_data

    @contextlib.contextmanager
    def test_context(self, runnables, test_timeout, keep_alive=False):
        self.events.append(("test_context_enter", keep_alive))
        try:
            yield self.test_data
        finally:
            self.events.append("test_context_exit")

    def serve_test_results(self):
        self.events.append("serve_test_results")


def test_serve_is_rejected_before_an_unsupported_engine_runs():
    engine = _UnsupportedEngine(create_test_context())

    with (
        patch("planemo.engine.test.engine_context", return_value=contextlib.nullcontext(engine)),
        pytest.raises(click.UsageError, match="managed Galaxy engine"),
    ):
        run_tests(create_test_context(), [], engine="cwltool", serve=True)


def test_reports_are_written_before_the_managed_server_waits():
    events = []
    test_data = SimpleNamespace(structured_data={"version": "0.1", "tests": []})
    engine = _RecordingEngine(events, test_data)

    def handle_reports(*args, **kwds):
        events.append("handle_reports")
        return 7

    with (
        patch("planemo.engine.test.engine_context", return_value=contextlib.nullcontext(engine)),
        patch("planemo.engine.test.handle_reports_and_summary", side_effect=handle_reports),
    ):
        exit_code = run_tests(create_test_context(), [], engine="galaxy", serve=True, test_timeout=31)

    assert exit_code == 7
    assert events == [
        ("test_context_enter", True),
        "handle_reports",
        "serve_test_results",
        "test_context_exit",
    ]
