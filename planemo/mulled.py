"""Planemo specific utilities for dealing with mulled containers.

The extend Galaxy/galaxy-lib's features with planemo specific idioms.
"""
from __future__ import absolute_import

import os

from galaxy.tools.deps.mulled.mulled_build import (
    DEFAULT_CHANNELS,
    ensure_installed,
    InvolucroContext,
)

from planemo.io import shell


def build_involucro_context(ctx, **kwds):
    """Build a galaxy-lib CondaContext tailored to planemo use.

    Using planemo's common command-line/global config options.
    """
    involucro_path_default = os.path.join(ctx.workspace, "involucro")
    involucro_path = kwds.get("involucro_path", involucro_path_default)
    use_planemo_shell = kwds.get("use_planemo_shell_exec", True)
    shell_exec = shell if use_planemo_shell else None
    involucro_context = InvolucroContext(involucro_bin=involucro_path,
                                         shell_exec=shell_exec)
    if not ensure_installed(involucro_context, True):
        raise Exception("Failed to install involucro for Planemo.")
    return involucro_context


def build_mull_target_kwds(ctx, **kwds):
    """Adapt Planemo's CLI and workspace configuration to galaxy-lib's mulled_build options."""
    involucro_context = build_involucro_context(ctx, **kwds)
    channels = kwds.get("conda_ensure_channels", ",".join(DEFAULT_CHANNELS))

    return {
        'involucro_context': involucro_context,
        'channels': channels.split(","),
    }

__all__ = [
    "build_involucro_context",
]
