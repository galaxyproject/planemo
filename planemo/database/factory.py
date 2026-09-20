"""Create a DatabaseSource from supplied planemo configuration."""

import contextlib
from typing import Optional

from .interface import DatabaseSource
from .postgres import LocalPostgresDatabaseSource
from .postgres_docker import DockerPostgresDatabaseSource
from .postgres_singularity import SingularityPostgresDatabaseSource

DATABASE_SOURCE_CLASSES = {
    "postgres": LocalPostgresDatabaseSource,
    "postgres_docker": DockerPostgresDatabaseSource,
    "postgres_singularity": SingularityPostgresDatabaseSource,
}


def database_source_class(database_type):
    """Return the database source class registered for ``database_type``."""
    try:
        return DATABASE_SOURCE_CLASSES[database_type]
    except KeyError:
        raise Exception("Unknown database type [%s]." % database_type) from None


def is_managed_database_type(database_type):
    """Return whether ``database_type`` names a registered database source."""
    return database_type in DATABASE_SOURCE_CLASSES


def create_database_source(
    profile_directory: Optional[str] = None, for_database_commands: bool = False, **kwds
) -> DatabaseSource:
    """Return a :class:`planemo.database.interface.DatabaseSource` for configuration."""
    database_type = kwds.get("database_type", "auto")
    if database_type == "auto":
        raise Exception(
            "Managing a database server requires naming one - pass --database_type with "
            "postgres, postgres_docker or postgres_singularity."
        )

    source_class = database_source_class(database_type)
    source_class.validate_configuration(
        profile_directory=profile_directory,
        for_database_commands=for_database_commands,
        **kwds,
    )
    return source_class(profile_directory=profile_directory, **kwds)


def started_database_source(
    profile_directory: Optional[str] = None, for_database_commands: bool = False, **kwds
) -> DatabaseSource:
    """Construct and start a :class:`planemo.database.interface.DatabaseSource`."""
    database_source = create_database_source(
        profile_directory=profile_directory,
        for_database_commands=for_database_commands,
        **kwds,
    )
    database_source.start()
    return database_source


@contextlib.contextmanager
def database_source_context(profile_directory: Optional[str] = None, for_database_commands: bool = False, **kwds):
    """Yield a started database source and stop it when doing so preserves its data."""
    database_source = started_database_source(
        profile_directory=profile_directory,
        for_database_commands=for_database_commands,
        **kwds,
    )
    try:
        yield database_source
    finally:
        if not database_source.keep_running_after_database_commands:
            database_source.stop()


__all__ = (
    "create_database_source",
    "database_source_class",
    "database_source_context",
    "is_managed_database_type",
    "started_database_source",
)
