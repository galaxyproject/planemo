
``job_config_init`` command
========================================

This section is auto-generated from the help text for the planemo command
``job_config_init``. This help message can be generated with ``planemo job_config_init
--help``.

**Usage**::

    planemo job_config_init [OPTIONS] TOOL_PATH

**Help**

Initialize an small Galaxy job config file for hosted workflow runs.
**Options**::


      --open                          Open the file in your default editor after
                                      creation.
      --docker / --no_docker          Run Galaxy tools in Docker if enabled.
      --docker_cmd TEXT               Command used to launch docker (defaults to
                                      docker).
      --docker_sudo / --no_docker_sudo
                                      Flag to use sudo when running docker.
      --docker_host TEXT              Docker host to target when executing docker
                                      commands (defaults to localhost).
      --docker_sudo_cmd TEXT          sudo command to use when --docker_sudo is
                                      enabled (defaults to sudo).
      --docker_run_extra_arguments TEXT
                                      Extra arguments to pass to docker run.
      --singularity / --no_singularity
                                      Run Galaxy tools in Singularity if enabled.
      --singularity_cmd TEXT          Command used to execute singularity (defaults
                                      to 'singularity').
      --singularity_sudo / --no_singularity_sudo
                                      Flag to use sudo when running docker.
      --singularity_sudo_cmd TEXT     sudo command to use when --singularity_sudo is
                                      enabled (defaults to sudo).
      --extra_tools PATH              Extra tool sources to include in Galaxy's tool
                                      panel (file or directory). These will not be
                                      linted/tested/etc... but they will be
                                      available to workflows and for interactive
                                      use.
      --test_data DIRECTORY           test-data directory to for specified tool(s).
      --tpv / --no_tpv                Include TPV (Total Perspective Vortex)
                                      configuration and shared usegalaxy* database
                                      of tool cores and memory for allocation
                                      purposes.
      --runner [local|slurm|drmaa|k8s|condor]
                                      Galaxy runner (e.g. DRM) to target.
      --galaxy_version TEXT           Version of Galaxy to target for configuration
                                      (default 24.2).
      --help                          Show this message and exit.
    
