
``serve`` command
======================================

This section is auto-generated from the help text for the planemo command
``serve``. This help message can be generated with ``planemo serve
--help``.

**Usage**::

    planemo serve [OPTIONS] TOOL_PATH

**Help**

Launch a Galaxy instance with the specified tool in the tool panel.

The Galaxy tool panel will include just the referenced tool or tools (by
default all the tools in the current working directory) and the upload
tool.

planemo will search parent directories to see if any is a Galaxy instance
- but one can pick the Galaxy instance to use with the ``--galaxy_root``
option or force planemo to download a disposable instance with the
``--install_galaxy`` flag.

``planemo`` uses temporarily generated config files and environment
variables to attempt to shield this execution of Galaxy from manually
launched runs against that same Galaxy root - but this may not be bullet
proof yet so please careful and do not try this against production Galaxy
instances.

**Options**::


      --galaxy_root DIRECTORY         Root of development galaxy directory to
                                      execute command with.
      --port INTEGER                  Port to serve Galaxy on (default is 9090).
      --install_galaxy                Download and configure a disposable copy of
                                      Galaxy from github.
      --no_cache_galaxy               Skip caching of Galaxy source and
                                      dependencies obtained with --install_galaxy.
                                      Not caching this results in faster downloads
                                      (no git) - so is better on throw away
                                      instances such with TravisCI.
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
    
