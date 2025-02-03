"""Utilities to help linting various targets."""

import os
from typing import (
    Any,
    Dict,
    TYPE_CHECKING,
)

from galaxy.tool_util.lint import (
    LintContext,
    Linter,
)

from planemo.io import error
from planemo.xml import validation

if TYPE_CHECKING:
    from planemo.cli import PlanemoCliContext


def build_lint_args(ctx: "PlanemoCliContext", **kwds) -> Dict[str, Any]:
    """Handle common report, error, and skip linting arguments."""
    report_level = kwds.get("report_level", "all")
    fail_level = kwds.get("fail_level", "warn")
    skip = kwds.get("skip", ctx.global_config.get("lint_skip"))
    if skip is None:
        skip = []
    if isinstance(skip, list):
        skip_types = skip
    else:
        skip_types = [s.strip() for s in skip.split(",")]

    for skip_file in kwds.get("skip_file", []):
        with open(skip_file) as f:
            for line in f.readlines():
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                skip_types.append(line)

    linters = Linter.list_linters()
    invalid_skip_types = list(set(skip_types) - set(linters))
    if len(invalid_skip_types):
        error(f"Unknown linter type(s) {invalid_skip_types} in list of linters to be skipped. Known linters {linters}")

    lint_args: Dict[str, Any] = dict(
        level=report_level,
        fail_level=fail_level,
        skip_types=skip_types,
    )
    return lint_args


def setup_lint(ctx, **kwds):
    """Prepare lint_args and lint_ctx to begin linting a target."""
    lint_args = kwds.get("lint_args", None) or build_lint_args(ctx, **kwds)
    lint_ctx = LintContext(level=lint_args["level"], skip_types=lint_args["skip_types"])
    return lint_args, lint_ctx


def handle_lint_complete(lint_ctx, lint_args, failed=False):
    """Complete linting of a target and decide exit code."""
    if not failed:
        failed = lint_ctx.failed(lint_args["fail_level"])
    if failed:
        error("Failed linting")
    return 1 if failed else 0


def lint_xsd(lint_ctx, schema_path, path):
    """Lint XML at specified path with supplied schema."""
    name = lint_ctx.object_name or os.path.basename(path)
    validator = validation.get_validator(require=True)
    validation_result = validator.validate(schema_path, path)
    if not validation_result.passed:
        msg = "Invalid XML found in file: %s. Errors [%s]"
        msg = msg % (name, validation_result.output)
        lint_ctx.error(msg)
    else:
        lint_ctx.info("File validates against XML schema.")


__all__ = (
    "build_lint_args",
    "handle_lint_complete",
    "lint_xsd",
)
