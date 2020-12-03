
``conda_install`` command
======================================

This section is auto-generated from the help text for the planemo command
``conda_install``. This help message can be generated with ``planemo conda_install
--help``.

**Usage**::

    planemo conda_install [OPTIONS] TARGET

**Help**

Install conda packages for tool requirements.
**Options**::


      -r, --recursive                 Recursively perform command for
                                      subdirectories.
    
      --conda_prefix DIRECTORY        Conda prefix to use for conda dependency
                                      commands.
    
      --conda_exec FILE               Location of conda executable.
      --conda_debug                   Enable more verbose conda logging.
      --conda_channels, --conda_ensure_channels TEXT
                                      Ensure conda is configured with specified
                                      comma separated list of channels.
    
      --conda_use_local               Use locally built packages while building
                                      Conda environments.
    
      --global                        Install Conda dependencies globally instead of
                                      in requirement specific environments packaged
                                      for tools. If the Conda bin directory is on
                                      your PATH, tools may still use binaries but
                                      this is more designed for interactive testing
                                      and debugging.
    
      --conda_auto_init / --no_conda_auto_init
                                      Conda dependency resolution for Galaxy will
                                      auto install conda itself using miniconda if
                                      not availabe on conda_prefix.
    
      --help                          Show this message and exit.
    
