
``conda_init`` command
======================================

This section is auto-generated from the help text for the planemo command
``conda_init``. This help message can be generated with ``planemo conda_init
--help``.

**Usage**::

    planemo conda_init [OPTIONS]

**Help**

Download and install conda.

This will download conda for managing dependencies for your platform
using the appropriate Miniconda installer.

By running this command, you are agreeing to the terms of the conda
license a 3-clause BSD 3 license. Please review full license at
http://docs.continuum.io/anaconda/eula.

Planemo will print a warning and terminate with an exit code of 7
if Conda is already installed.

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
    
