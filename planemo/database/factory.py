"""Create a DatabaseSource from supplied planemo configuration."""
from .postgres import LocalPostgresDatabaseSource


def create_database_source(**kwds):
    """Return a :class:`planemo.database.DatabaseSource` for configuration."""
    database_type = kwds.get("database_type", "postgres")
    if database_type == "postgres":
        return LocalPostgresDatabaseSource(**kwds)
    # TODO
    # from .sqlite import SqliteDatabaseSource
    # elif database_type == "sqlite":
    #     return SqliteDatabaseSource(**kwds)
    else:
        raise Exception("Unknown database type [%s]." % database_type)


__all__ = (
    "create_database_source",
)
