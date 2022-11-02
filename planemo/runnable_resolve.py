import os

import requests

from planemo.galaxy.profiles import translate_alias
from planemo.galaxy.workflows import GALAXY_WORKFLOWS_PREFIX
from planemo.tools import uri_to_path
from .runnable import (
    for_path,
    for_uri,
    GALAXY_TOOLS_PREFIX,
)


def for_runnable_identifier(ctx, runnable_identifier, kwds):
    """Convert URI, path, or alias into Runnable."""
    # could be a URI, path, or alias
    current_profile = kwds.get("profile")
    runnable_identifier = translate_alias(ctx, runnable_identifier, current_profile)
    if not runnable_identifier.startswith(GALAXY_WORKFLOWS_PREFIX):
        runnable_identifier = uri_to_path(ctx, runnable_identifier)
    if os.path.exists(runnable_identifier):
        runnable = for_path(runnable_identifier)
    else:  # assume galaxy workflow or tool id
        if "/repos/" in runnable_identifier:
            runnable_identifier = f"{GALAXY_TOOLS_PREFIX}{runnable_identifier}"
        elif not runnable_identifier.startswith(GALAXY_WORKFLOWS_PREFIX):
            runnable_identifier = f"{GALAXY_WORKFLOWS_PREFIX}{runnable_identifier}"
        runnable = for_uri(runnable_identifier)
    return runnable


def for_runnable_identifiers(ctx, runnable_identifiers, kwds):
    """Convert lists of URIs, paths, and/or aliases into Runnables."""
    runnables = []
    for r in runnable_identifiers:
        runnable = for_runnable_identifier(ctx, r, kwds)
        if isinstance(runnable, list):
            runnables.extend(runnable)
        else:
            runnables.append(runnable)
    return runnables


def install_args_list_to_runnables(ctx, install_args_list, kwds):
    runnables = []
    for repo in install_args_list:
        url = f'{repo["tool_shed_url"]}api/repositories/get_repository_revision_install_info'
        response = requests.get(
            url, params={"name": repo["name"], "owner": repo["owner"], "changeset_revision": repo["changeset_revision"]}
        )
        response.raise_for_status()
        install_info = response.json()
        repository_metadata = install_info[1]
        assert repository_metadata["model_class"] == "RepositoryMetadata"
        for tool in repository_metadata.get("valid_tools", []):
            runnable = for_runnable_identifier(ctx, tool["guid"], kwds)
            runnables.append(runnable)
    return runnables
