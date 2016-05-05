
``tool_factory`` command
======================================

This section is auto-generated from the help text for the planemo command
``tool_factory``. This help message can be generated with ``planemo tool_factory
--help``.

**Usage**::

    planemo tool_factory [OPTIONS]

**Help**

(Experimental) Launch Galaxy with the Tool Factory 2 available.

For more information about the Galaxy Tool Factory see the publication
Creating reusable tools from scripts: the Galaxy Tool Factory by Lazarus
et. al. (10.1093/bioinformatics/bts573). Available at
http://www.ncbi.nlm.nih.gov/pubmed/23024011.

**Options**::


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
      --extra_tools PATH              Extra tool sources to include in Galaxy's tool
                                      panel (file or directory). These will not be
                                      linted/tested/etc... but they will be
                                      available to workflows and for interactive
                                      use.
      --daemon                        Serve Galaxy process as a daemon.
      --pid_file TEXT                 Location of pid file is executed with
                                      --daemon.
      --help                          Show this message and exit.
    
