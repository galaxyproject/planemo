"""Manage a PostgreSQL server in a Singularity/Apptainer container."""

import os
import subprocess
import time
from tempfile import mkdtemp
from typing import Optional

from galaxy.util.commands import shell_process

from planemo.io import (
    info,
    terminate_process_group,
)
from .interface import (
    DatabaseConfigurationError,
    DatabaseSource,
)
from .postgres import (
    _CommandBuilder,
    ExecutesPostgresSqlMixin,
)

DEFAULT_POSTGRES_DATABASE_NAME = "galaxy"
DEFAULT_POSTGRES_USER = "galaxy"
DEFAULT_POSTGRES_PASSWORD = "mysecretpassword"
DEFAULT_DOCKERIMAGE = "postgres:14.2-alpine3.15"
DEFAULT_STARTUP_TIMEOUT = 120
DEFAULT_STOP_TIMEOUT = 15
CONTAINER_SOCKET_DIRECTORY = "/var/run/postgresql"
POSTGRES_SOCKET_NAME = ".s.PGSQL.5432"


def start_postgres_singularity(
    singularity_command,
    database_location,
    databasename=DEFAULT_POSTGRES_DATABASE_NAME,
    user=DEFAULT_POSTGRES_USER,
    password=DEFAULT_POSTGRES_PASSWORD,
):
    info(f"Postgres database stored at: {database_location}")
    pgdata_path = os.path.join(database_location, "pgdata")
    pgrun_path = os.path.join(database_location, "pgrun")
    log_path = os.path.join(database_location, "postgres.log")

    if not os.path.exists(pgdata_path):
        os.makedirs(pgdata_path)
    if not os.path.exists(pgrun_path):
        os.makedirs(pgrun_path)

    run_command = [
        *singularity_command,
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
        "POSTGRES_INITDB_ARGS=--encoding=UTF-8",
        f"docker://{DEFAULT_DOCKERIMAGE}",
    ]
    info("Starting postgres singularity container")
    with open(log_path, "ab", buffering=0) as log:
        return shell_process(
            run_command,
            start_new_session=True,
            stdout=log,
            stderr=subprocess.STDOUT,
        )


class SingularityPostgresDatabaseSource(ExecutesPostgresSqlMixin, DatabaseSource):
    """Postgres database source running inside a Singularity container."""

    store_connection_in_profile = False
    PROFILE_OPTIONS = (
        "postgres_storage_location",
        "singularity_cmd",
        "singularity_sudo",
        "singularity_sudo_cmd",
    )

    def __init__(self, profile_directory: Optional[str] = None, **kwds):
        """Construct a postgres database source from planemo configuration."""

        singularity_cmd = kwds.get("singularity_cmd") or "singularity"
        self.singularity_command = []
        if kwds.get("singularity_sudo", False):
            self.singularity_command.append(kwds.get("singularity_sudo_cmd") or "sudo")
        self.singularity_command.append(singularity_cmd)
        self.database_user = DEFAULT_POSTGRES_USER
        self.database_password = DEFAULT_POSTGRES_PASSWORD
        if kwds.get("postgres_storage_location") is not None:
            self.database_location = kwds["postgres_storage_location"]
        elif profile_directory:
            self.database_location = os.path.join(profile_directory, "postgres")
        else:
            self.database_location = mkdtemp(suffix="_planemo_postgres_db")
        self.database_location = os.path.abspath(os.path.expanduser(self.database_location))
        self.database_socket_dir = os.path.join(self.database_location, "pgrun")
        self.log_file = os.path.join(self.database_location, "postgres.log")
        self.startup_timeout = DEFAULT_STARTUP_TIMEOUT
        self.stop_timeout = DEFAULT_STOP_TIMEOUT
        self._kwds = kwds
        self.running_process = None

    def start(self):
        if self.running_process is not None and self.running_process.poll() is None:
            return
        self.running_process = start_postgres_singularity(
            singularity_command=self.singularity_command,
            database_location=self.database_location,
            user=self.database_user,
            password=self.database_password,
        )
        deadline = time.monotonic() + self.startup_timeout
        try:
            while True:
                return_code = self.running_process.poll()
                if return_code is not None:
                    self.running_process = None
                    raise RuntimeError(
                        f"PostgreSQL Singularity container exited during startup with code {return_code}; "
                        f"see {self.log_file}."
                    )
                if self._database_is_ready():
                    return
                if time.monotonic() >= deadline:
                    raise RuntimeError(
                        f"PostgreSQL Singularity container did not become ready within {self.startup_timeout} seconds; "
                        f"see {self.log_file}."
                    )
                time.sleep(1)
        except BaseException:
            self.stop()
            raise

    def stop(self):
        process = self.running_process
        if process is None:
            return
        self.running_process = None
        if process.poll() is not None:
            process.wait()
            return
        if terminate_process_group(process.pid, timeout=self.stop_timeout, reap=process.poll):
            process.wait()
        else:
            info("PostgreSQL Singularity container process group could not be stopped.")

    def _singularity_exec_command(self, executable, *args):
        return [
            *self.singularity_command,
            "exec",
            "-B",
            f"{self.database_socket_dir}:{CONTAINER_SOCKET_DIRECTORY}",
            "-e",
            "-C",
            f"docker://{DEFAULT_DOCKERIMAGE}",
            executable,
            *args,
        ]

    def _database_is_ready(self):
        socket_path = os.path.join(self.database_socket_dir, POSTGRES_SOCKET_NAME)
        if not os.path.exists(socket_path):
            return False
        command = self._singularity_exec_command(
            "pg_isready",
            "--host",
            CONTAINER_SOCKET_DIRECTORY,
            "--username",
            self.database_user,
            "--dbname",
            "postgres",
        )
        try:
            completed = subprocess.run(
                command,
                check=False,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                timeout=10,
            )
        except subprocess.TimeoutExpired:
            return False
        return completed.returncode == 0

    def _psql_command_builder(self, *args):
        command_builder = _CommandBuilder(*self._singularity_exec_command("psql"))
        command_builder.append_command("--tuples-only")
        command_builder.append_command("--username", self.database_user)
        command_builder.append_command("--host", CONTAINER_SOCKET_DIRECTORY)
        command_builder.append_command("--dbname", "postgres")
        command_builder.append_command("-P", "pager=off")
        command_builder.extend_command(args)
        return command_builder

    def sqlalchemy_url(self, identifier):
        """Return URL for PostgreSQL connection via Unix socket."""
        return "postgresql://%s:%s@/%s?host=%s" % (
            self.database_user,
            self.database_password,
            identifier,
            self.database_socket_dir,
        )

    @classmethod
    def validate_configuration(cls, profile_directory=None, for_database_commands=False, **kwds):
        """Require stable storage for separate database administration commands."""
        if for_database_commands and profile_directory is None and not kwds.get("postgres_storage_location"):
            raise DatabaseConfigurationError(
                "Database administration with postgres_singularity requires --postgres-storage-location so separate "
                "commands operate on the same PostgreSQL cluster."
            )

    def profile_options(self):
        """Return enough configuration to restart this profile's container."""
        options = super().profile_options()
        options["postgres_storage_location"] = self.database_location
        return options


__all__ = ("SingularityPostgresDatabaseSource",)
