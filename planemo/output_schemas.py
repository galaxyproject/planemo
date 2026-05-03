"""JSON Schema exports for Planemo machine-readable outputs."""

from __future__ import annotations

from typing import Any

from planemo import __version__
from planemo.output_models import (
    PlanemoInvocationDownloadManifest,
    PlanemoRunOutputs,
    PlanemoTestReport,
)

SCHEMA_VERSION = "0.1"

SCHEMA_MODELS = {
    "invocation-download-manifest": PlanemoInvocationDownloadManifest,
    "run-outputs": PlanemoRunOutputs,
    "test-report": PlanemoTestReport,
}


def load_output_schemas(schema_name: str | None = None) -> dict[str, Any]:
    """Return Planemo output JSON Schemas."""
    schemas = {}
    for name, model in SCHEMA_MODELS.items():
        if schema_name is None or name == schema_name:
            schema = model.model_json_schema()
            schema["$schema"] = "https://json-schema.org/draft/2020-12/schema"
            schemas[name] = schema
    if schema_name is not None and schema_name not in schemas:
        raise KeyError(schema_name)
    return {
        "schema_version": SCHEMA_VERSION,
        "planemo_version": __version__,
        "schemas": schemas,
    }


__all__ = (
    "SCHEMA_MODELS",
    "SCHEMA_VERSION",
    "load_output_schemas",
)
