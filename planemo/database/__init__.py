"""Package is responsible for managing planemo profile databases.

This package makes it very easy to create and destroy databases, therefore it
should not be used for production data - and should not even be connnected
to a production database server.
"""

from .factory import (
    create_database_source,
    database_source_class,
    database_source_context,
    is_managed_database_type,
    started_database_source,
)
from .interface import DatabaseConfigurationError

__all__ = (
    "create_database_source",
    "database_source_class",
    "database_source_context",
    "DatabaseConfigurationError",
    "is_managed_database_type",
    "started_database_source",
)
