
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
      --galaxy_sqlite_database DIRECTORY
                                      Preseeded Galaxy sqlite database to target.
      --install_galaxy                Download and configure a disposable copy of
                                      Galaxy from github.
      --no_cache_galaxy               Skip caching of Galaxy source and dependencies
                                      obtained with --install_galaxy. Not caching
                                      this results in faster downloads (no git) - so
                                      is better on throw away instances such with
                                      TravisCI.
      --no_cleanup                    Do not cleanup temp files created for and by
                                      Galaxy.
      --job_config_file PATH          Job configuration file for Galaxy to target.
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
      --tool_dependency_dir DIRECTORY
                                      Tool dependency dir for Galaxy to target.
      --brew_dependency_resolution    Configure Galaxy to use plain brew dependency
                                      resolution.
      --shed_dependency_resolution    Configure Galaxy to use brewed Tool Shed
                                      dependency resolution.
      --help                          Show this message and exit.
    
