
``profile_delete`` command
======================================

This section is auto-generated from the help text for the planemo command
``profile_delete``. This help message can be generated with ``planemo profile_delete
--help``.

**Usage**::

    planemo profile_delete [OPTIONS] PROFILE_NAME

**Help**

Delete a profile.
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
    
      --docker_cmd TEXT               Command used to launch docker (defaults to
                                      docker).
    
      --docker_sudo / --no_docker_sudo
                                      Flag to use sudo when running docker.
      --docker_host TEXT              Docker host to target when executing docker
                                      commands (defaults to localhost).
    
      --docker_sudo_cmd TEXT          sudo command to use when --docker_sudo is
                                      enabled (defaults to sudo).
    
      --help                          Show this message and exit.
    
