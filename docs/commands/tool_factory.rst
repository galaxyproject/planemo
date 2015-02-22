
``tool_factory`` command
===============================

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
      --install_galaxy                Download and configure a disposable copy of
                                      Galaxy from github.
      --no_cleanup                    Do not cleanup temp files created for and by
                                      Galaxy.
      --test_data DIRECTORY           test-data directory to for specified
                                      tool(s).
      --dependency_resolvers_config_file PATH
                                      Dependency resolver configuration for Galaxy
                                      to target.
      --job_config_file PATH          Job configuration file for Galaxy to target.
      --tool_dependency_dir DIRECTORY
                                      Tool dependency dir for Galaxy to target.
      --brew_dependency_resolution    Configure Galaxy to use plain brew
                                      dependency resolution.
      --help                          Show this message and exit.
    
