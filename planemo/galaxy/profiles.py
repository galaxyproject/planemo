"""This modules describes the abstraction of a Galaxy profile.

This is a workspace with a specific default configuration and shed
tool setup. It is meant to be used with various serve commands.
"""
import json
import os
import shutil

from .config import (
    DATABASE_LOCATION_TEMPLATE,
    attempt_database_preseed,
)
from planemo.database import create_database_source

PROFILE_OPTIONS_JSON_NAME = "planemo_profile_options.json"
ALREADY_EXISTS_EXCEPTION = "Cannot create profile with name [%s], directory [%s] already exists."


def profile_exists(ctx, profile_name, **kwds):
    """Return a truthy value iff the specified profile already exists."""
    profile_directory = _profile_directory(ctx, profile_name)
    return os.path.exists(profile_directory)


def list_profiles(ctx, **kwds):
    return os.listdir(ctx.galaxy_profiles_directory)


def delete_profile(ctx, profile_name, **kwds):
    """Delete profile with the specified name."""
    profile_directory = _profile_directory(ctx, profile_name)
    profile_options = _read_profile_options(profile_directory)
    database_type = profile_options.get("database_type")
    if database_type != "sqlite":
        database_source = create_database_source(**kwds)
        database_identifier = _profile_to_database_identifier(profile_name)
        database_source.delete_database(
            database_identifier,
        )
    shutil.rmtree(profile_directory)


def create_profile(ctx, profile_name, **kwds):
    """Create a profile with the specified name."""
    profile_directory = _profile_directory(ctx, profile_name)
    if profile_exists(ctx, profile_name, **kwds):
        message = ALREADY_EXISTS_EXCEPTION % (
            profile_name, profile_directory
        )
        raise Exception(message)

    os.makedirs(profile_directory)
    database_type = kwds.get("database_type", "sqlite")
    if database_type != "sqlite":
        database_source = create_database_source(**kwds)
        database_identifier = _profile_to_database_identifier(profile_name)
        database_source.create_database(
            database_identifier,
        )
        database_connection = database_source.sqlalchemy_url(database_identifier)
    else:
        database_location = os.path.join(profile_directory, "galaxy.sqlite")
        attempt_database_preseed(None, database_location, **kwds)
        database_connection = DATABASE_LOCATION_TEMPLATE % database_location

    stored_profile_options = {
        'database_type': database_type,
        'database_connection': database_connection,
    }
    profile_options_path = _stored_profile_options_path(profile_directory)
    with open(profile_options_path, "w") as f:
        json.dump(stored_profile_options, f)


def ensure_profile(ctx, profile_name, **kwds):
    """Ensure a Galaxy profile exists and return profile defaults."""
    if not profile_exists(ctx, profile_name, **kwds):
        create_profile(ctx, profile_name, **kwds)

    return _profile_options(ctx, profile_name, **kwds)


def _profile_options(ctx, profile_name, **kwds):
    profile_directory = _profile_directory(ctx, profile_name)
    file_path = os.path.join(profile_directory, "files")
    shed_tool_path = os.path.join(profile_directory, "shed_tools")
    shed_tool_conf = os.path.join(profile_directory, "shed_tool_conf.xml")
    tool_dependency_dir = os.path.join(profile_directory, "deps")

    profile_options = _read_profile_options(profile_directory)
    profile_options.update(dict(
        file_path=file_path,
        tool_dependency_dir=tool_dependency_dir,
        shed_tool_conf=shed_tool_conf,
        shed_tool_path=shed_tool_path,
    ))

    return profile_options


def _profile_to_database_identifier(profile_name):
    char_lst = [c if c.isalnum() else "_" for c in profile_name]
    return "plnmoprof_%s" % "".join(char_lst)


def _read_profile_options(profile_directory):
    profile_options_path = _stored_profile_options_path(profile_directory)
    with open(profile_options_path, "r") as f:
        profile_options = json.load(f)
    return profile_options


def _stored_profile_options_path(profile_directory):
    profile_options_path = os.path.join(
        profile_directory, PROFILE_OPTIONS_JSON_NAME
    )
    return profile_options_path


def _profile_directory(ctx, profile_name):
    return os.path.join(ctx.galaxy_profiles_directory, profile_name)

__all__ = [
    "ensure_profile",
    "create_profile",
    "delete_profile",
]
