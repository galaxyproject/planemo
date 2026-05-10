"""Pydantic models for Planemo machine-readable test outputs."""

from __future__ import annotations

from typing import (
    Any,
    Literal,
)

from pydantic import (
    BaseModel,
    ConfigDict,
    field_validator,
    model_validator,
)

PlanemoTestStatus = Literal["success", "failure", "error", "skip"]


class PlanemoTestSummary(BaseModel):
    model_config = ConfigDict(extra="allow")

    num_tests: int
    num_failures: int
    num_skips: int
    num_errors: int


class PlanemoTestCaseData(BaseModel):
    model_config = ConfigDict(extra="allow")

    status: PlanemoTestStatus
    inputs: dict[str, Any] | None = None
    job: dict[str, Any] | None = None
    invocation_details: dict[str, Any] | None = None
    problem_log: str | None = None
    output_problems: list[str] | None = None
    execution_problem: str | None = None
    start_datetime: str | None = None
    end_datetime: str | None = None

    @field_validator("status", mode="before")
    @classmethod
    def normalize_legacy_status(cls, value):
        if value == "skipped":
            return "skip"
        return value


class PlanemoTestCase(BaseModel):
    model_config = ConfigDict(extra="allow")

    id: str
    has_data: bool
    data: PlanemoTestCaseData | None = None
    doc: str | None = None
    test_type: str | None = None

    @model_validator(mode="after")
    def require_data_when_present(self):
        if self.has_data and self.data is None:
            raise ValueError("data is required when has_data is true")
        return self


class PlanemoTestReport(BaseModel):
    model_config = ConfigDict(extra="allow")

    version: str = "0.1"
    tests: list[PlanemoTestCase]
    summary: PlanemoTestSummary | None = None
    exit_code: int | None = None


__all__ = (
    "PlanemoTestCase",
    "PlanemoTestCaseData",
    "PlanemoTestReport",
    "PlanemoTestStatus",
    "PlanemoTestSummary",
)
