
``shed_test`` command
========================================

This section is auto-generated from the help text for the planemo command
``shed_test``. This help message can be generated with ``planemo shed_test
--help``.

**Usage**::

    planemo shed_test [OPTIONS] PROJECT

**Help**

Run tests of published shed artifacts.

This command will start a Galaxy instance configured to target the
specified shed, find published artifacts (tools and dependencies)
corresponding to command-line arguments and ``.shed.yml`` file(s),
install these artifacts, and run the tool tests for these commands.

**Options**::


      -r, --recursive                 Recursively perform command for nested
                                      repository directories.
      --shed_fail_fast                If multiple repositories are specified and an
                                      error occurs stop immediately instead of
                                      processing remaining repositories.
      --owner TEXT                    Tool Shed repository owner (username).
      --name TEXT                     Tool Shed repository name (defaults to the
                                      inferred tool directory name).
      --shed_email TEXT               E-mail for Tool Shed auth (required unless
                                      shed_key is specified).
      --shed_key TEXT                 API key for Tool Shed access. An API key is
                                      required unless e-mail and password is
                                      specified. This key can be specified with
                                      either --shed_key or --shed_key_from_env.
      --shed_key_from_env TEXT        Environment variable to read API key for Tool
                                      Shed access from.
      --shed_password TEXT            Password for Tool Shed auth (required unless
                                      shed_key is specified).
      -t, --shed_target TEXT          Tool Shed to target (this can be 'toolshed',
                                      'testtoolshed', 'local' (alias for
                                      http://localhost:9009/), an arbitrary url or
                                      mappings defined ~/.planemo.yml.
      --galaxy_root DIRECTORY         Root of development galaxy directory to
                                      execute command with.
      --galaxy_python_version [3|3.8|3.9|3.10|3.11|3.12]
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
      --docker_run_extra_arguments TEXT
                                      Extra arguments to pass to docker run.
      --mulled_containers, --biocontainers
                                      Test tools against mulled containers (forces
                                      --docker). Disables conda resolution unless
                                      any conda option has been set explicitly.
      --galaxy_startup_timeout INTEGER RANGE
                                      Wait for galaxy to start before assuming
                                      Galaxy did not start.  [x>=1]
      --job_config_file FILE          Job configuration file for Galaxy to target.
      --tool_dependency_dir DIRECTORY
                                      Tool dependency dir for Galaxy to target.
      --tool_data_path DIRECTORY      Directory where data used by tools is located.
                                      Required if tests are run in docker and should
                                      make use of external reference data.
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
      --conda_channels, --conda_ensure_channels TEXT
                                      Ensure conda is configured with specified
                                      comma separated list of channels.
      --conda_use_local               Use locally built packages while building
                                      Conda environments.
      --conda_dependency_resolution   Configure Galaxy to use only conda for
                                      dependency resolution.
      --conda_auto_install / --no_conda_auto_install
                                      Conda dependency resolution for Galaxy will
                                      attempt to install requested but missing
                                      packages.
      --conda_auto_init / --no_conda_auto_init
                                      Conda dependency resolution for Galaxy will
                                      auto install conda itself using miniforge if
                                      not availabe on conda_prefix.
      --simultaneous_uploads / --no_simultaneous_uploads
                                      When uploading files to Galaxy for tool or
                                      workflow tests or runs, upload multiple files
                                      simultaneously without waiting for the
                                      previous file upload to complete.
      --check_uploads_ok / --no_check_uploads_ok
                                      When uploading files to Galaxy for tool or
                                      workflow tests or runs, check that the history
                                      is in an 'ok' state before beginning tool or
                                      workflow execution.
      --profile TEXT                  Name of profile (created with the
                                      profile_create command) to use with this
                                      command.
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
      --file_path DIRECTORY           Location for files created by Galaxy (e.g.
                                      database/files).
      --database_connection TEXT      Database connection string to use for Galaxy.
      --postgres-storage-location TEXT
                                      storage path for postgres database, used for
                                      local singularity postgres.
      --shed_tool_conf TEXT           Location of shed tools conf file for Galaxy.
      --shed_tool_path TEXT           Location of shed tools directory for Galaxy.
      --galaxy_single_user / --no_galaxy_single_user
                                      By default Planemo will configure Galaxy to
                                      run in single-user mode where there is just
                                      one user and this user is automatically logged
                                      it. Use --no_galaxy_single_user to prevent
                                      Galaxy from running this way.
      --paste_test_data_paths / --no_paste_test_data_paths
                                      By default Planemo will use or not use
                                      Galaxy's path paste option to load test data
                                      into a history based on the engine type it is
                                      targeting. This can override the logic to
                                      explicitly enable or disable path pasting.
      --update_test_data              Update test-data directory with job outputs
                                      (normally written to directory
                                      --job_output_files if specified.)
      --test_output PATH              Output test report (HTML - for humans)
                                      defaults to tool_test_output.html.
      --test_output_text PATH         Output test report (Basic text - for display
                                      in CI)
      --test_output_markdown PATH     Output test report (Markdown style - for
                                      humans & computers)
      --test_output_markdown_minimal PATH
                                      Output test report (Minimal markdown style -
                                      jost the table)
      --test_output_xunit PATH        Output test report (xunit style - for CI
                                      systems
      --test_output_junit PATH        Output test report (jUnit style - for CI
                                      systems
      --test_output_allure DIRECTORY  Output test allure2 framework resutls
      --test_output_json PATH         Output test report (planemo json) defaults to
                                      tool_test_output.json.
      --job_output_files DIRECTORY    Write job outputs to specified directory.
      --summary [none|minimal|compact]
                                      Summary style printed to planemo's standard
                                      output (see output reports for more complete
                                      summary). Set to 'none' to disable completely.
      --test_timeout INTEGER          Maximum runtime of a single test in seconds.
      --fail_fast                     Stop on first job failure.
      --skip_dependencies             Do not install shed dependencies as part of
                                      repository installation.
      --help                          Show this message and exit.
    
