
``database_create`` command
========================================

This section is auto-generated from the help text for the planemo command
``database_create``. This help message can be generated with ``planemo database_create
--help``.

**Usage**::

    planemo database_create [OPTIONS] IDENTIFIER

**Help**

Create a *development* database.

Currently the only implementation is postgres which will be managed with
``psql``.

Planemo ``database_`` commands make it very easy to create and destroy
databases, therefore it should not be used for production data - and it
should not even be connnected to a production database server. Planemo
is intended for development purposes only.

Planemo will assume that it can manage and access postgres databases
without specifying a password. This can be accomplished by configuring
postgres to not required a password for the planemo user or by specifying
a password in a ``.pgpass`` file.

Planemo can be configured to not require a password for the planemo user in
the postgres configuration file ``pg_hba.conf`` (on Ubuntu/Debian linux
distros this file is in /etc/postgresql/<postgres_version>/main/ directory).
Adding the following lines to that file will allow planemo and Galaxy to
access the databases without a password.


# "local" is for Unix domain socket connections only
local   all   all                    trust
# IPv4 local connections:
host    all   all    127.0.0.1/32    trust
# IPv6 local connections:
host    all   all    ::1/128         trust

More information on the ``pg_hda.conf`` configuration file can be found at
http://www.postgresql.org/docs/9.3/static/auth-pg-hba-conf.html.

Information on ``.pgpass`` files can be found at at the following location:
http://www.postgresql.org/docs/9.4/static/libpq-pgpass.html. In Ubuntu and
Debian distros - a postgres user likely already exists and its password can
be set by setting up a file ``~/.pgpass`` file with the following contents.


*:*:*:postgres:<postgres_password>

**Options**::


      --postgres                      Use postgres database type.
      --database_type [postgres|postgres_docker|postgres_singularity|sqlite|auto]
                                      Type of database to use for profile - 'auto',
                                      'sqlite', 'postgres', 'postgres_docker' , and
                                      postgres_singularity are available options.
                                      Use postgres to use an existing postgres
                                      server you user can access without a password
                                      via the psql command. Use postgres_docker to
                                      have Planemo manage a docker container running
                                      postgres. . Use  postgres_singularity to have
                                      Planemo run postgres using
                                      singularity/apptainer. Data with
                                      postgres_docker is not yet persisted past when
                                      you restart the docker container launched by
                                      Planemo so be careful with this option.
      --postgres_psql_path TEXT       Name or or path to postgres client binary
                                      (psql).
      --postgres_database_user TEXT   Postgres username for managed development
                                      databases.
      --postgres_database_host TEXT   Postgres host name for managed development
                                      databases.
      --postgres_database_port TEXT   Postgres port for managed development
                                      databases.
      --docker_cmd TEXT               Command used to launch docker (defaults to
                                      docker).
      --docker_sudo / --no_docker_sudo
                                      Flag to use sudo when running docker.
      --docker_host TEXT              Docker host to target when executing docker
                                      commands (defaults to localhost).
      --docker_sudo_cmd TEXT          sudo command to use when --docker_sudo is
                                      enabled (defaults to sudo).
      --docker_run_extra_arguments TEXT
                                      Extra arguments to pass to docker run.
      --help                          Show this message and exit.
    
