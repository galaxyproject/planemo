import os

from planemo.galaxy.profiles import translate_alias
from planemo.galaxy.workflows import GALAXY_WORKFLOWS_PREFIX
from planemo.tools import uri_to_path
from .runnable import (
    for_path,
    for_uri,
)


def for_runnable_identifier(ctx, runnable_identifier, kwds, temp_path=None, return_all=False):
    """Convert URI, path, or alias into Runnable."""
    # could be a URI, path, or alias
    current_profile = kwds.get("profile")
    runnable_identifier = translate_alias(ctx, runnable_identifier, current_profile)
    if not runnable_identifier.startswith(GALAXY_WORKFLOWS_PREFIX):
        runnable_identifier = uri_to_path(ctx, runnable_identifier)
    if os.path.exists(runnable_identifier):
        runnable = for_path(runnable_identifier, temp_path=temp_path, return_all=return_all)
    else:  # assume galaxy workflow id
        if not runnable_identifier.startswith(GALAXY_WORKFLOWS_PREFIX):
            runnable_identifier = f"{GALAXY_WORKFLOWS_PREFIX}{runnable_identifier}"
        runnable = for_uri(runnable_identifier)
    return runnable


def for_runnable_identifiers(ctx, runnable_identifiers, kwds, temp_path=None):
    """Convert lists of URIs, paths, and/or aliases into Runnables."""
    runnables = []
    for r in runnable_identifiers:
        runnable = for_runnable_identifier(ctx, r, kwds, temp_path=temp_path, return_all=True)
        if isinstance(runnable, list):
            runnables.extend(runnable)
        else:
            runnables.append(runnable)
    return runnables
