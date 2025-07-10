"""Module describes a :class:`DatabaseSource` for managed, dockerized postgres databases."""

import os
import shutil
import time
from tempfile import mkdtemp
from typing import Optional

from galaxy.util.commands import shell_process

from planemo.io import info
from .interface import DatabaseSource
from .postgres import ExecutesPostgresSqlMixin

DEFAULT_CONTAINER_NAME = "planemopostgres"
DEFAULT_POSTGRES_DATABASE_NAME = "galaxy"
DEFAULT_POSTGRES_USER = "galaxy"
DEFAULT_POSTGRES_PASSWORD = "mysecretpassword"
DEFAULT_DOCKERIMAGE = "postgres:14.2-alpine3.15"


def start_postgres_singularity(
    singularity_path,
    database_location,
    databasename=DEFAULT_POSTGRES_DATABASE_NAME,
    user=DEFAULT_POSTGRES_USER,
    password=DEFAULT_POSTGRES_PASSWORD,
):
    info(f"Postgres database stored at: {database_location}")
    pgdata_path = os.path.join(database_location, "pgdata")
    pgrun_path = os.path.join(database_location, "pgrun")

    if not os.path.exists(pgdata_path):
        os.makedirs(pgdata_path)
    if not os.path.exists(pgrun_path):
        os.makedirs(pgrun_path)

    VERSION_FILE = os.path.join(pgdata_path, "PG_VERSION")
    run_command = [
        singularity_path,
        "run",
        "-B",
        f"{pgdata_path}:/var/lib/postgresql/data",
        "-B",
        f"{pgrun_path}:/var/run/postgresql",
        "-e",
        "-C",
        "--env",
        f"POSTGRES_DB={databasename}",
        "--env",
        f"POSTGRES_USER={user}",
        "--env",
        f"POSTGRES_PASSWORD={password}",
        "--env",
        "POSTGRES_INITDB_ARGS='--encoding=UTF-8'",
        f"docker://{DEFAULT_DOCKERIMAGE}",
    ]
    info("Starting postgres singularity container")
    p = shell_process(run_command)
    # Give the container time to initialize the database
    for _ in range(10):
        if os.path.exists(VERSION_FILE):
            break
        time.sleep(5)
        info("Waiting for the postgres database to initialize.")
    else:
        try:
            p.terminate()
        except Exception as e:
            info(f"Failed to terminate process: {e}")
        raise Exception("Failed to initialize the postgres database.")
    return p


class SingularityPostgresDatabaseSource(ExecutesPostgresSqlMixin, DatabaseSource):
    """
    Postgres database running inside a Singularity container. Should be used with
    "with" statements to automatically start and stop the container.
    """

    def __init__(self, profile_directory: Optional[str] = None, **kwds):
        """Construct a postgres database source from planemo configuration."""

        self.singularity_path = "singularity"
        self.database_user = DEFAULT_POSTGRES_USER
        self.database_password = DEFAULT_POSTGRES_PASSWORD
        if kwds.get("postgres_storage_location") is not None:
            self.database_location = kwds["postgres_storage_location"]
        elif profile_directory:
            self.database_location = os.path.join(profile_directory, "postgres")
        else:
            self.database_location = os.path.join(mkdtemp(suffix="_planemo_postgres_db"))
        self.database_socket_dir = os.path.join(self.database_location, "pgrun")
        self.container_instance_name = f"{DEFAULT_CONTAINER_NAME}-{int(time.time() * 1000000)}"
        self._kwds = kwds
        self.running_process = None

    def create_database(self, identifier):
        # Not needed, we'll create the database automatically when the container starts.
        pass

    def delete_database(self, identifier):
        shutil.rmtree(self.database_location, ignore_errors=True)

    def __enter__(self):
        self.running_process = start_postgres_singularity(
            singularity_path=self.singularity_path,
            database_location=self.database_location,
            user=self.database_user,
            password=self.database_password,
        )
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        if self.running_process:
            try:
                self.running_process.terminate()
                postmaster_pid_file = os.path.join(self.database_location, "pgdata", "postmaster.pid")
                if os.path.exists(postmaster_pid_file):
                    os.remove(postmaster_pid_file)
                postmaster_lock_file = os.path.join(self.database_location, "pgrun", "pgrun/.s.PGSQL.5432.lock")
                if os.path.exists(postmaster_lock_file):
                    os.remove(postmaster_lock_file)
            except Exception as e:
                info(f"Failed to terminate process: {e}")

    def sqlalchemy_url(self, identifier):
        """Return URL for PostgreSQL connection via Unix socket."""
        return "postgresql://%s:%s@/%s?host=%s" % (
            self.database_user,
            self.database_password,
            identifier,
            self.database_socket_dir,
        )


__all__ = ("SingularityPostgresDatabaseSource",)
