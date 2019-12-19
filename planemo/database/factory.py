"""Create a DatabaseSource from supplied planemo configuration."""
from galaxy.tool_util.deps.commands import which

from .postgres import LocalPostgresDatabaseSource
from .postgres_docker import DockerPostgresDatabaseSource


def create_database_source(**kwds):
    """Return a :class:`planemo.database.DatabaseSource` for configuration."""
    database_type = kwds.get("database_type", "auto")
    if database_type == "auto":
        if which("psql"):
            database_type = "postgres"
        elif which("docker"):
            database_type = "postgres_docker"
        else:
            raise Exception("Cannot find executables for psql or docker, cannot configure a database source.")

    if database_type == "postgres":
        return LocalPostgresDatabaseSource(**kwds)
    elif database_type == "postgres_docker":
        return DockerPostgresDatabaseSource(**kwds)
    # TODO
    # from .sqlite import SqliteDatabaseSource
    # elif database_type == "sqlite":
    #     return SqliteDatabaseSource(**kwds)
    else:
        raise Exception("Unknown database type [%s]." % database_type)


__all__ = (
    "create_database_source",
)
