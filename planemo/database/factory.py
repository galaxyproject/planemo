"""Create a DatabaseSource from supplied planemo configuration."""
from .postgres import LocalPostgresDatabaseSource


def create_database_source(**kwds):
    """Return a :class:`planemo.database.DatabaseSource` for configuration."""
    return LocalPostgresDatabaseSource(**kwds)


__all__ = ['create_database_source']
