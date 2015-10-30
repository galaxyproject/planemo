""" Tool linting module that lints Galaxy tools for their URLs
"""

import planemo.lint


def lint_tool_urls(root, lint_ctx):
    planemo.lint.lint_urls(root, lint_ctx)
