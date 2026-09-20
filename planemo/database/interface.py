"""Describe the interface classes of the planemo.database package."""

import abc
from typing import Optional


class DatabaseConfigurationError(ValueError):
    """Indicate that a database source cannot support the requested operation."""


class DatabaseSource(metaclass=abc.ABCMeta):
    """Interface describing a source of profile databases."""

    keep_running_after_database_commands = False
    store_connection_in_profile = True
    PROFILE_OPTIONS: tuple[str, ...] = ()

    @abc.abstractmethod
    def create_database(self, identifier):
        """Create a database with specified short identifier.

        Throw an exception if it already exists.
        """

    @abc.abstractmethod
    def delete_database(self, identifier):
        """Delete a database with specified short identifier.

        Throw an exception if it does not exist.
        """

    @abc.abstractmethod
    def list_databases(self):
        """Return identifiers associated with database source."""

    @abc.abstractmethod
    def sqlalchemy_url(self, identifier) -> Optional[str]:
        """Return a URL string for use by sqlalchemy."""

    def start(self):
        """Start the database source, if necessary."""
        pass

    def stop(self):
        """Stop the database source, if necessary."""
        pass

    @classmethod
    def validate_configuration(cls, profile_directory=None, for_database_commands=False, **kwds):
        """Validate backend-specific configuration before allocating resources."""

    def profile_options(self):
        """Return backend configuration that a persistent profile must retain."""
        kwds = getattr(self, "_kwds", {})
        return {option: kwds[option] for option in self.PROFILE_OPTIONS if option in kwds}


__all__ = (
    "DatabaseConfigurationError",
    "DatabaseSource",
)
