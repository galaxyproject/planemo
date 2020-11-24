
``conda_env`` command
======================================

This section is auto-generated from the help text for the planemo command
``conda_env``. This help message can be generated with ``planemo conda_env
--help``.

**Usage**::

    planemo conda_env [OPTIONS] TOOL_PATH

**Help**

Activate a conda environment for tool.

Source the output of this command to activate a conda environment for this
tool.

::

    $ . <(planemo conda_env seqtk_seq.xml)
    Deactivate environment with conda_env_deactivate
    (seqtk_seq_v6) $ which seqtk
    /home/planemo/miniconda2/envs/jobdepsDkzcjjfecc6d406196737781ff4456ec60975c137e04884e4f4b05dc68192f7cec4656/bin/seqtk
    (seqtk_seq_v6) $ conda_env_deactivate
    $


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
    
