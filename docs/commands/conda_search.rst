
``conda_search`` command
======================================

This section is auto-generated from the help text for the planemo command
``conda_search``. This help message can be generated with ``planemo conda_search
--help``.

**Usage**::

    planemo conda_search [OPTIONS] TERM

**Help**

Perform conda search with Planemo's conda.

Implicitly adds channels Planemo is configured with.

**Options**::


      --conda_prefix DIRECTORY        Conda prefix to use for conda dependency
                                      commands.
    
      --conda_exec FILE               Location of conda executable.
      --conda_debug                   Enable more verbose conda logging.
      --conda_channels, --conda_ensure_channels TEXT
                                      Ensure conda is configured with specified
                                      comma separated list of channels.
    
      --conda_use_local               Use locally built packages while building
                                      Conda environments.
    
      --help                          Show this message and exit.
    
