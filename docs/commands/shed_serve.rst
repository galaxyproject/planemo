
``shed_serve`` command
======================================

This section is auto-generated from the help text for the planemo command
``shed_serve``. This help message can be generated with ``planemo shed_serve
--help``.

**Usage**::

    planemo shed_serve [OPTIONS] PROJECT

**Help**

Serve a transient Galaxy with published repositories installed.

This command will start a Galaxy instance configured to target the
specified shed, find published artifacts (tools and dependencies)
corresponding to command-line arguments and ``.shed.yml`` file(s),
install these artifacts, and serve a Galaxy instances that can be
logged into and explored interactively.

**Options**::


      -r, --recursive                 Recursively perform command for nested
                                      repository directories.
      --fail_fast                     If multiple repositories are specified and an
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
      --galaxy_database_seed PATH     Preseeded Galaxy sqlite database to target.
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
      --job_config_file PATH          Job configuration file for Galaxy to target.
      --tool_dependency_dir DIRECTORY
                                      Tool dependency dir for Galaxy to target.
      --port INTEGER                  Port to serve Galaxy on (default is 9090).
      --host TEXT                     Host to bind Galaxy to. Default is 127.0.0.1
                                      that is restricted to localhost connections
                                      for security reasons set to 0.0.0.0 to bind
                                      Galaxy to all ports including potentially
                                      publicly accessible ones.
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
      --conda_prefix DIRECTORY        Conda prefix to use for conda dependency
                                      commands.
      --conda_exec PATH               Location of conda executable.
      --conda_debug                   Enable more verbose conda logging.
      --conda_ensure_channels TEXT    Ensure conda is configured with specified
                                      comma separated list of channels.
      --conda_dependency_resolution   Configure Galaxy to use only conda for
                                      dependency resolution.
      --conda_copy_dependencies       Conda dependency resolution for Galaxy will
                                      copy dependencies instead of attempting to
                                      link them.
      --conda_auto_install            Conda dependency resolution for Galaxy will
                                      auto install will attempt to install requested
                                      but missing packages.
      --conda_auto_init               Conda dependency resolution for Galaxy will
                                      auto install conda itself using miniconda if
                                      not availabe on conda_prefix.
      --profile TEXT                  Location of pid file is executed with
                                      --daemon.
      --file_path DIRECTORY           Location for files created by Galaxy (e.g.
                                      database/files).
      --database_connection TEXT      Database connection string to use for Galaxy.
      --shed_tool_conf TEXT           Location of shed tools conf file for Galaxy.
      --shed_tool_path TEXT           Location of shed tools directory for Galaxy.
      --pid_file TEXT                 Location of pid file is executed with
                                      --daemon.
      --skip_dependencies             Do not install shed dependencies as part of
                                      repository installation.
      --help                          Show this message and exit.
    
