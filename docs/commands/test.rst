
``test`` command
======================================

This section is auto-generated from the help text for the planemo command
``test``. This help message can be generated with ``planemo test
--help``.

**Usage**::

    planemo test [OPTIONS] TOOL_PATH

**Help**

Run specified tool's tests within Galaxy.

All referenced tools (by default all the tools in the current working
directory) will be tested and the results quickly summarized.

To run these tests planemo needs a Galaxy instance to utilize, planemo
will search parent directories to see if any is a Galaxy instance
- but one can pick the Galaxy instance to use with the --galaxy_root
option or force planemo to download a disposable instance with the
``--install_galaxy`` flag.

In additon to to quick summary printed to the console - various detailed
output summaries can be configured. ``tool_test_output.html`` (settable
via ``--test_output``) will contain a human consumable HTML report
describing the test run. A JSON file (settable via ``--test_output_json``
and defaulting to ``tool_test_output.json``) will also be created. These
files can can be disabled by passing in empty arguments or globally by
setting the values ``default_test_output`` and/or
``default_test_output_json`` in ``~/.planemo.yml`` to ``null``. For
continuous integration testing a xUnit-style report can be confiured using
the ``--test_output_xunit``.

planemo uses temporarily generated config files and environment variables
to attempt to shield this execution of Galaxy from manually launched runs
against that same Galaxy root - but this may not be bullet proof yet so
please careful and do not try this against production Galaxy instances.

**Options**::


      --failed                        Re-run only failed tests. This command will
                                      read tool_test_output.json to determine which
                                      tests failed so this file must have been
                                      produced with the same set of tool ids
                                      previously.
      --galaxy_root DIRECTORY         Root of development galaxy directory to
                                      execute command with.
      --galaxy_database_seed PATH     Preseeded Galaxy sqlite database to target.
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
                                      environment for Galaxy, this should be used or
                                      instance to preserve an externally configured
                                      virtual environment or conda environment.
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
      --job_config_file PATH          Job configuration file for Galaxy to target.
      --tool_dependency_dir DIRECTORY
                                      Tool dependency dir for Galaxy to target.
      --test_data DIRECTORY           test-data directory to for specified tool(s).
      --tool_data_table PATH          tool_data_table_conf.xml file to for specified
                                      tool(s).
      --dependency_resolvers_config_file PATH
                                      Dependency resolver configuration for Galaxy
                                      to target.
      --brew_dependency_resolution    Configure Galaxy to use plain brew dependency
                                      resolution.
      --shed_dependency_resolution    Configure Galaxy to use brewed Tool Shed
                                      dependency resolution.
      --no_dependency_resolution      Configure Galaxy with no dependency resolvers.
      --conda_prefix DIRECTORY        Conda prefix to use for conda dependency
                                      commands.
      --conda_exec PATH               Location of conda executable.
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
      --update_test_data              Update test-data directory with job outputs
                                      (normally written to directory
                                      --job_output_files if specified.)
      --test_output PATH              Output test report (HTML - for humans)
                                      defaults to tool_test_output.html.
      --test_output_text PATH         Output test report (Basic text - for display
                                      in CI)
      --test_output_markdown PATH     Output test report (Markdown style - for
                                      humans & computers)
      --test_output_xunit PATH        Output test report (xUnit style - for
                                      computers).
      --test_output_json PATH         Output test report (planemo json) defaults to
                                      tool_test_output.json.
      --job_output_files DIRECTORY    Write job outputs to specified directory.
      --summary [none|minimal|compact]
                                      Summary style printed to planemo's standard
                                      output (see output reports for more complete
                                      summary). Set to 'none' to disable completely.
      --engine [galaxy|docker_galaxy|cwltool]
                                      Select an engine to run or test aritfacts such
                                      as tools and workflows. Defaults to a local
                                      Galaxy, but running Galaxy within a Docker
                                      container or the CWL reference implementation
                                      'cwltool' and be selected.
      --non_strict_cwl                Disable strict validation of CWL.
      --no-container, --no_container  If cwltool engine is used, disable Docker
                                      container usage.
      --docker_galaxy_image TEXT      Docker image identifier for docker-galaxy-
                                      flavor used if engine type is specified as
                                      ``docker-galaxy``. Defaults to to bgruening
                                      /galaxy-stable.
      --help                          Show this message and exit.
    
