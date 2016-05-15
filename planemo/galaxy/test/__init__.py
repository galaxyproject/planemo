"""Entry point and interface for ``planemo.galaxy.test`` package."""
from .actions import run_in_config
from .structures import StructuredData
from .actions import handle_reports
from .actions import handle_reports_and_summary


__all__ = [
    "handle_reports",
    "handle_reports_and_summary",
    "run_in_config",
    "StructuredData",
]
