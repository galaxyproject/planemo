
``training_generate_from_wf`` command
======================================

This section is auto-generated from the help text for the planemo command
``training_generate_from_wf``. This help message can be generated with ``planemo training_generate_from_wf
--help``.

**Usage**::

    planemo training_generate_from_wf [OPTIONS] TOOL_PATH

**Help**

Create tutorial skeleton from workflow.
**Options**::


      --topic_name TEXT               Name (directory name) of the topic to create
                                      or in which a tutorial should be created or
                                      updates  [required]
    
      --tutorial_name TEXT            Name (directory name) of the tutorial to
                                      modify  [required]
    
      --workflow PATH                 Workflow of the tutorial (locally)
      --galaxy_url TEXT               URL of a Galaxy instance with the workflow
      --galaxy_api_key TEXT           API key on the Galaxy instance with the
                                      workflow
    
      --workflow_id TEXT              ID of the workflow on the Galaxy instance
      --galaxy_root DIRECTORY         Root of development galaxy directory to
                                      execute command with.
    
      --galaxy_python_version [3|3.6|3.7|3.8|3.9]
                                      Python version to start Galaxy under
      --extra_tools PATH              Extra tool sources to include in Galaxy's tool
                                      panel (file or directory). These will not be
                                      linted/tested/etc... but they will be
                                      available to workflows and for interactive
                                      use.
    
      --install_galaxy                Download and configure a disposable copy of
                                      Galaxy from github.
    
      --galaxy_branch TEXT            Branch of Galaxy to target (defaults to
                                      master) if a Galaxy root isn't specified.
    
      --galaxy_source TEXT            Git source of Galaxy to target (defaults to
                                      the official galaxyproject github source if a
                                      Galaxy root isn't specified.
    
      --skip_venv                     Do not create or source a virtualenv
                                      environment for Galaxy, this should be used to
                                      preserve an externally configured virtual
                                      environment or conda environment.
    
      --no_cache_galaxy               Skip caching of Galaxy source and dependencies
                                      obtained with --install_galaxy. Not caching
                                      this results in faster downloads (no git) - so
                                      is better on throw away instances such with
                                      TravisCI.
    
      --no_cleanup                    Do not cleanup temp files created for and by
                                      Galaxy.
    
      --galaxy_email TEXT             E-mail address to use when launching single-
                                      user Galaxy server.
    
      --docker / --no_docker          Run Galaxy tools in Docker if enabled.
      --docker_cmd TEXT               Command used to launch docker (defaults to
                                      docker).
    
      --docker_sudo / --no_docker_sudo
                                      Flag to use sudo when running docker.
      --docker_host TEXT              Docker host to target when executing docker
                                      commands (defaults to localhost).
    
      --docker_sudo_cmd TEXT          sudo command to use when --docker_sudo is
                                      enabled (defaults to sudo).
    
      --mulled_containers, --biocontainers
                                      Test tools against mulled containers (forces
                                      --docker).
    
      --job_config_file FILE          Job configuration file for Galaxy to target.
      --tool_dependency_dir DIRECTORY
                                      Tool dependency dir for Galaxy to target.
      --port INTEGER                  Port to serve Galaxy on (default is 9090).
      --host TEXT                     Host to bind Galaxy to. Default is 127.0.0.1
                                      that is restricted to localhost connections
                                      for security reasons set to 0.0.0.0 to bind
                                      Galaxy to all ports including potentially
                                      publicly accessible ones.
    
      --engine [galaxy|docker_galaxy|external_galaxy]
                                      Select an engine to serve artifacts such as
                                      tools and workflows. Defaults to a local
                                      Galaxy, but running Galaxy within a Docker
                                      container.
    
      --non_strict_cwl                Disable strict validation of CWL.
      --docker_galaxy_image TEXT      Docker image identifier for docker-galaxy-
                                      flavor used if engine type is specified as
                                      ``docker-galaxy``. Defaults to
                                      quay.io/bgruening/galaxy.
    
      --docker_extra_volume PATH      Extra path to mount if --engine docker.
      --test_data DIRECTORY           test-data directory to for specified tool(s).
      --tool_data_table PATH          tool_data_table_conf.xml file to for specified
                                      tool(s).
    
      --dependency_resolvers_config_file FILE
                                      Dependency resolver configuration for Galaxy
                                      to target.
    
      --brew_dependency_resolution    Configure Galaxy to use plain brew dependency
                                      resolution.
    
      --shed_dependency_resolution    Configure Galaxy to use brewed Tool Shed
                                      dependency resolution.
    
      --no_dependency_resolution      Configure Galaxy with no dependency resolvers.
      --conda_prefix DIRECTORY        Conda prefix to use for conda dependency
                                      commands.
    
      --conda_exec FILE               Location of conda executable.
      --conda_debug                   Enable more verbose conda logging.
      --conda_channels, --conda_ensure_channels TEXT
                                      Ensure conda is configured with specified
                                      comma separated list of channels.
    
      --conda_use_local               Use locally built packages while building
                                      Conda environments.
    
      --conda_dependency_resolution   Configure Galaxy to use only conda for
                                      dependency resolution.
    
      --conda_copy_dependencies       Conda dependency resolution for Galaxy will
                                      copy dependencies instead of attempting to
                                      link them.
    
      --conda_auto_install / --no_conda_auto_install
                                      Conda dependency resolution for Galaxy will
                                      attempt to install requested but missing
                                      packages.
    
      --conda_auto_init / --no_conda_auto_init
                                      Conda dependency resolution for Galaxy will
                                      auto install conda itself using miniconda if
                                      not availabe on conda_prefix.
    
      --profile TEXT                  Name of profile (created with the
                                      profile_create command) to use with this
                                      command.
    
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
    
      --file_path DIRECTORY           Location for files created by Galaxy (e.g.
                                      database/files).
    
      --database_connection TEXT      Database connection string to use for Galaxy.
      --shed_tool_conf TEXT           Location of shed tools conf file for Galaxy.
      --shed_tool_path TEXT           Location of shed tools directory for Galaxy.
      --galaxy_single_user / --no_galaxy_single_user
                                      By default Planemo will configure Galaxy to
                                      run in single-user mode where there is just
                                      one user and this user is automatically logged
                                      it. Use --no_galaxy_single_user to prevent
                                      Galaxy from running this way.
    
      --daemon                        Serve Galaxy process as a daemon.
      --pid_file FILE                 Location of pid file is executed with
                                      --daemon.
    
      --ignore_dependency_problems    When installing shed repositories for
                                      workflows, ignore dependency issues. These
                                      likely indicate a problem but in some cases
                                      may not prevent a workflow from successfully
                                      executing.
    
      --skip_client_build             Do not build Galaxy client when serving
                                      Galaxy.
    
      --shed_install / --no_shed_install
                                      By default Planemo will attempt to install
                                      repositories needed for workflow testing. This
                                      may not be appropriate for production servers
                                      and so this can disabled by calling planemo
                                      with --no_shed_install.
    
      --help                          Show this message and exit.
    
