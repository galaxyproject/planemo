
``conda_env`` command
======================================

This section is auto-generated from the help text for the planemo command
``conda_env``. This help message can be generated with ``planemo conda_env
--help``.

**Usage**::

    planemo conda_env [OPTIONS] TOOL_PATH

**Help**

Source output to activate a conda environment for this tool.

    % . <(planemo conda_env bowtie2.xml)
    % which bowtie2
    TODO_PLACE_PATH_HERE

**Options**::


      --conda_prefix DIRECTORY      Conda prefix to use for conda dependency
                                    commands.
      --conda_exec PATH             Location of conda executable.
      --conda_debug                 Enable more verbose conda logging.
      --conda_ensure_channels TEXT  Ensure conda is configured with specified comma
                                    separated list of channels.
      --help                        Show this message and exit.
    
