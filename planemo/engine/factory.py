"""Module contains factory method for building class:`Engine` objects."""

import contextlib

from .galaxy import GalaxyEngine
from .cwltool import CwlToolEngine

UNKNOWN_ENGINE_TYPE_MESSAGE = "Unknown engine type specified [%s]."


def is_galaxy_engine(**kwds):
    """Return True iff the engine configured is :class:`GalaxyEngine`."""
    engine_type_str = kwds.get("engine", "galaxy")
    return engine_type_str == "galaxy"


def build_engine(ctx, **kwds):
    """Build an engine from the supplied planemo configuration."""
    engine_type_str = kwds.get("engine", "galaxy")

    if engine_type_str == "galaxy":
        engine_type = GalaxyEngine
    elif engine_type_str == "cwltool":
        engine_type = CwlToolEngine
    else:
        raise Exception(UNKNOWN_ENGINE_TYPE_MESSAGE % engine_type_str)

    return engine_type(ctx, **kwds)


@contextlib.contextmanager
def engine_context(ctx, **kwds):
    """A :func:`contextlib.contextmanager` engine builder for use with ``with`` statements.

    https://docs.python.org/2/library/contextlib.html
    """
    engine = None
    try:
        engine = build_engine(ctx, **kwds)
        yield engine
    finally:
        if engine is not None:
            engine.cleanup()


__all__ = [
    "is_galaxy_engine",
    "build_engine",
    "engine_context",
]
