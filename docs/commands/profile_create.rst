
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
      --database_type [postgres|postgres_docker|sqlite|auto]
                                      Type of database to use for profile - 'auto',
                                      'sqlite', 'postgres', and 'postgres_docker'
                                      are available options. Use postgres to use an
                                      existing postgres server you user can access
                                      without a password via the psql command. Use
                                      postgres_docker to have Planemo manage a
                                      docker container running postgres. Data with
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
    
      --engine [galaxy|docker_galaxy|external_galaxy]
                                      Select an engine to serve artifacts such as
                                      tools and workflows. Defaults to a local
                                      Galaxy, but running Galaxy within a Docker
                                      container.
    
      --docker_cmd TEXT               Command used to launch docker (defaults to
                                      docker).
    
      --docker_sudo / --no_docker_sudo
                                      Flag to use sudo when running docker.
      --docker_host TEXT              Docker host to target when executing docker
                                      commands (defaults to localhost).
    
      --docker_sudo_cmd TEXT          sudo command to use when --docker_sudo is
                                      enabled (defaults to sudo).
    
      --galaxy_url TEXT               Remote Galaxy URL to use with external Galaxy
                                      engine.
    
      --galaxy_user_key TEXT          User key to use with external Galaxy engine.
      --galaxy_admin_key TEXT         Admin key to use with external Galaxy engine.
      --help                          Show this message and exit.
    
