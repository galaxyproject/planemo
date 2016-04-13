"""This modules describes the abstraction of a Galaxy profile.

This is a workspace with a specific default configuration and shed
tool setup. It is meant to be used with various serve commands.
"""
import os

from .config import (
    DATABASE_LOCATION_TEMPLATE,
    attempt_database_preseed,
)


def ensure_profile(ctx, profile_name, **kwds):
    """Ensure a Galaxy profile exists and return profile defaults."""
    profile_directory = _profile_directory(ctx, profile_name)
    database_location = os.path.join(profile_directory, "galaxy.sqlite")

    if not os.path.exists(profile_directory):
        os.makedirs(profile_directory)
        attempt_database_preseed(None, database_location, **kwds)

    database_connection = DATABASE_LOCATION_TEMPLATE % database_location
    file_path = os.path.join(profile_directory, "files")
    shed_tool_path = os.path.join(profile_directory, "shed_tools")
    shed_tool_conf = os.path.join(profile_directory, "shed_tool_conf.xml")
    tool_dependency_dir = os.path.join(profile_directory, "deps")

    return dict(
        file_path=file_path,
        database_connection=database_connection,
        tool_dependency_dir=tool_dependency_dir,
        shed_tool_conf=shed_tool_conf,
        shed_tool_path=shed_tool_path,
    )


def _profile_directory(ctx, profile_name):
    return os.path.join(ctx.galaxy_profiles_directory, profile_name)

__all__ = [
    "ensure_profile",
]
