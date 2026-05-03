"""Pydantic models for Planemo machine-readable outputs."""

from __future__ import annotations

from typing import Any

from pydantic import RootModel

from planemo.test.models import PlanemoTestReport


class PlanemoRunOutputs(RootModel[dict[str, Any]]):
    """Permissive model for ``planemo run --output_json`` outputs."""


__all__ = (
    "PlanemoRunOutputs",
    "PlanemoTestReport",
)
