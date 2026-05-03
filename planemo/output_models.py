"""Pydantic models for Planemo machine-readable outputs."""

from __future__ import annotations

from typing import (
    Any,
    Literal,
)

from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    RootModel,
)

from planemo.test.models import PlanemoTestReport


class PlanemoRunOutputs(RootModel[dict[str, Any]]):
    """Permissive model for ``planemo run --output_json`` outputs."""


def validated_run_outputs(data):
    return PlanemoRunOutputs.model_validate(data).model_dump(mode="json")


def validated_test_report(data):
    return PlanemoTestReport.model_validate(data).model_dump(mode="json")


class PlanemoInvocationDownloadMissingOutput(BaseModel):
    model_config = ConfigDict(extra="allow")

    id: str
    reason: Literal["missing", "skipped"] = "missing"


class PlanemoInvocationDownloadManifest(BaseModel):
    model_config = ConfigDict(extra="allow")

    invocation_id: str
    output_directory: str
    path_type: Literal["relative", "absolute"] = "relative"
    outputs: dict[str, Any]
    missing_outputs: list[PlanemoInvocationDownloadMissingOutput] = Field(default_factory=list)
    output_json: str | None = None


__all__ = (
    "PlanemoInvocationDownloadManifest",
    "PlanemoInvocationDownloadMissingOutput",
    "PlanemoRunOutputs",
    "PlanemoTestReport",
    "validated_run_outputs",
    "validated_test_report",
)
