"""Unit tests for --failed/--lf (last-failed) test filtering."""

import json
import os
from tempfile import NamedTemporaryFile
from typing import (
    Any,
    Dict,
    List,
)
from unittest.mock import MagicMock

import pytest

import planemo.engine.interface as _engine_iface
from planemo.engine.interface import BaseEngine
from planemo.test.results import StructuredData

# Importing test_utils first causes planemo.cli → planemo.runnable to be fully
# initialized before planemo.engine.interface is imported below, preventing the
# circular import that otherwise occurs when planemo.engine is loaded in isolation.
from .test_utils import create_test_context

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_json(tests: List[Dict[str, Any]]) -> str:
    tmp = NamedTemporaryFile(suffix=".json", delete=False, mode="w")
    json.dump({"version": "0.1", "tests": tests}, tmp)
    tmp.close()
    return tmp.name


def _test_entry(test_id: str, status: str) -> Dict[str, Any]:
    return {"id": test_id, "has_data": True, "data": {"status": status}}


def _make_test_case(test_id: str, index: int) -> MagicMock:
    tc = MagicMock()
    tc._test_id = test_id
    tc.index = index
    return tc


class _MinimalEngine(BaseEngine):
    """Stub engine that records which test cases reach _collect_test_results."""

    handled_runnable_types = []

    def __init__(self, ctx, submitted, **kwds):
        super().__init__(ctx, **kwds)
        self._submitted = submitted

    def _check_can_run_all(self, runnables):
        """No-op: allow any runnable type in tests."""

    def _run(self, runnables, job_paths, output_collectors=None):
        return []

    def _collect_test_results(self, cases, test_timeout):
        self._submitted.extend(cases)
        return []


def _run_engine_filter(test_cases, kwds):
    """Call BaseEngine.test() with a stub cases() and return submitted cases."""
    submitted = []
    engine = _MinimalEngine(create_test_context(), submitted, **kwds)

    fake_runnable = MagicMock()
    original_cases = _engine_iface.cases
    _engine_iface.cases = lambda _: test_cases
    try:
        engine.test([fake_runnable], test_timeout=0)
    finally:
        _engine_iface.cases = original_cases

    return submitted


# ---------------------------------------------------------------------------
# StructuredData.failed_ids
# ---------------------------------------------------------------------------


class TestFailedIds:
    def test_returns_only_non_success(self):
        json_path = _make_json(
            [
                _test_entry("tool_a_0", "success"),
                _test_entry("tool_a_1", "failure"),
                _test_entry("tool_b_0", "error"),
                _test_entry("tool_c_0", "skip"),
            ]
        )
        try:
            assert StructuredData(json_path=json_path).failed_ids == {"tool_a_1", "tool_b_0", "tool_c_0"}
        finally:
            os.unlink(json_path)

    def test_empty_when_all_pass(self):
        json_path = _make_json(
            [
                _test_entry("tool_a_0", "success"),
                _test_entry("tool_a_1", "success"),
            ]
        )
        try:
            assert StructuredData(json_path=json_path).failed_ids == set()
        finally:
            os.unlink(json_path)

    def test_all_when_all_fail(self):
        json_path = _make_json(
            [
                _test_entry("tool_a_0", "failure"),
                _test_entry("tool_b_0", "failure"),
            ]
        )
        try:
            assert StructuredData(json_path=json_path).failed_ids == {"tool_a_0", "tool_b_0"}
        finally:
            os.unlink(json_path)


# ---------------------------------------------------------------------------
# BaseEngine.test() filtering via --failed
# ---------------------------------------------------------------------------


class TestBaseEngineFailedFilter:
    def test_no_filter_without_flag(self, tmp_path):
        """Without --failed, all test cases are submitted."""
        cases = [_make_test_case("tool_a", 0), _make_test_case("tool_a", 1)]
        submitted = _run_engine_filter(cases, {"failed": False})
        assert len(submitted) == 2

    def test_filters_to_failed_only(self, tmp_path):
        """With --failed, only test cases matching failed IDs are submitted."""
        json_path = str(tmp_path / "prev.json")
        with open(json_path, "w") as f:
            json.dump(
                {
                    "version": "0.1",
                    "tests": [
                        _test_entry("tool_a_0", "failure"),
                        _test_entry("tool_a_1", "success"),
                    ],
                },
                f,
            )

        cases = [_make_test_case("tool_a", 0), _make_test_case("tool_a", 1)]
        submitted = _run_engine_filter(cases, {"failed": True, "test_output_json": json_path})
        assert len(submitted) == 1
        assert submitted[0]._test_id == "tool_a"
        assert submitted[0].index == 0

    def test_failed_json_takes_precedence_over_test_output_json(self, tmp_path):
        """--failed_json is used instead of --test_output_json when both are set."""
        output_json = str(tmp_path / "output.json")
        with open(output_json, "w") as f:
            json.dump(
                {
                    "version": "0.1",
                    "tests": [
                        _test_entry("tool_a_0", "success"),
                        _test_entry("tool_a_1", "failure"),
                    ],
                },
                f,
            )

        failed_json = str(tmp_path / "failed.json")
        with open(failed_json, "w") as f:
            json.dump(
                {
                    "version": "0.1",
                    "tests": [
                        _test_entry("tool_a_0", "failure"),
                        _test_entry("tool_a_1", "success"),
                    ],
                },
                f,
            )

        cases = [_make_test_case("tool_a", 0), _make_test_case("tool_a", 1)]
        submitted = _run_engine_filter(
            cases,
            {"failed": True, "failed_json": failed_json, "test_output_json": output_json},
        )
        # failed_json says index 0 failed; output.json says index 1 failed
        assert len(submitted) == 1
        assert submitted[0].index == 0

    def test_raises_when_no_json_exists(self, tmp_path):
        """--failed without a readable JSON raises ClickException."""
        from click import ClickException

        cases = [_make_test_case("tool_a", 0)]
        with pytest.raises(ClickException, match="--failed/--lf requires"):
            _run_engine_filter(cases, {"failed": True, "test_output_json": str(tmp_path / "missing.json")})

    def test_returns_empty_when_all_passed(self, tmp_path):
        """--failed when no tests previously failed returns zero submitted cases."""
        json_path = str(tmp_path / "prev.json")
        with open(json_path, "w") as f:
            json.dump(
                {
                    "version": "0.1",
                    "tests": [
                        _test_entry("tool_a_0", "success"),
                    ],
                },
                f,
            )

        cases = [_make_test_case("tool_a", 0)]
        submitted = _run_engine_filter(cases, {"failed": True, "test_output_json": json_path})
        assert submitted == []
