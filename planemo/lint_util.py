"""Utilities to help linting various targets."""
from __future__ import absolute_import

from galaxy.tools.lint import LintContext

from planemo.io import error
from planemo.tool_lint import (
    build_lint_args,  # TODO: Move to here.
)


def setup_lint(ctx, **kwds):
    """Setup lint_args and lint_ctx to begin linting a target."""
    lint_args = build_lint_args(ctx, **kwds)
    lint_ctx = LintContext(lint_args["level"])
    return lint_args, lint_ctx


def handle_lint_complete(lint_ctx, lint_args, failed=False):
    """Complete linting of a target and decide exit code."""
    if not failed:
        failed = lint_ctx.failed(lint_args["fail_level"])
    if failed:
        error("Failed linting")
    return 1 if failed else 0


__all__ = [
    "handle_lint_complete"
]
