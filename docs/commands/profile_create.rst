
``profile_create`` command
======================================

This section is auto-generated from the help text for the planemo command
``profile_create``. This help message can be generated with ``planemo profile_create
--help``.

**Usage**::

    planemo profile_create [OPTIONS] PROFILE_NAME

**Help**

Create a profile.
**Options**::


      --postgres                      Use postgres database type.
      --database_type [postgres|sqlite]
                                      Type of database to use for profile -
                                      currently only 'postgres' is available.
      --postgres_psql_path TEXT       Name or or path to postgres client binary
                                      (psql).
      --postgres_database_user TEXT   Postgres username for managed development
                                      databases.
      --postgres_database_host TEXT   Postgres host name for managed development
                                      databases.
      --postgres_database_port TEXT   Postgres port for managed development
                                      databases.
      --engine [galaxy|docker_galaxy]
                                      Select an engine to serve aritfacts such as
                                      tools and workflows. Defaults to a local
                                      Galaxy, but running Galaxy within a Docker
                                      container.
      --help                          Show this message and exit.
    
