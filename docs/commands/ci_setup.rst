
``ci_setup`` command
========================================

This section is auto-generated from the help text for the planemo command
``ci_setup``. This help message can be generated with ``planemo ci_setup
--help``.

**Usage**::

    planemo ci_setup [OPTIONS]

**Help**


Launch Galaxy instance, then terminate instance.

Useful for populating a CI cache.

**Options**::


      --galaxy_root DIRECTORY         Root of development galaxy directory to
                                      execute command with.
      --galaxy_python_version [3|3.7|3.8|3.9|3.10|3.11]
                                      Python version to start Galaxy under
      --extra_tools PATH              Extra tool sources to include in Galaxy's tool
                                      panel (file or directory). These will not be
                                      linted/tested/etc... but they will be
                                      available to workflows and for interactive
                                      use.
      --install_galaxy                Download and configure a disposable copy of
                                      Galaxy from github.
      --galaxy_branch TEXT            Branch of Galaxy to target (defaults to
                                      master) if a Galaxy root isn't specified.
      --galaxy_source TEXT            Git source of Galaxy to target (defaults to
                                      the official galaxyproject github source if a
                                      Galaxy root isn't specified.
      --skip_venv                     Do not create or source a virtualenv
                                      environment for Galaxy, this should be used to
                                      preserve an externally configured virtual
                                      environment or conda environment.
      --no_cache_galaxy               Skip caching of Galaxy source and dependencies
                                      obtained with --install_galaxy. Not caching
                                      this results in faster downloads (no git) - so
                                      is better on throw away instances such with
                                      TravisCI.
      --no_cleanup                    Do not cleanup temp files created for and by
                                      Galaxy.
      --galaxy_email TEXT             E-mail address to use when launching single-
                                      user Galaxy server.
      --docker / --no_docker          Run Galaxy tools in Docker if enabled.
      --docker_cmd TEXT               Command used to launch docker (defaults to
                                      docker).
      --docker_sudo / --no_docker_sudo
                                      Flag to use sudo when running docker.
      --docker_host TEXT              Docker host to target when executing docker
                                      commands (defaults to localhost).
      --docker_sudo_cmd TEXT          sudo command to use when --docker_sudo is
                                      enabled (defaults to sudo).
      --mulled_containers, --biocontainers
                                      Test tools against mulled containers (forces
                                      --docker). Disables conda resolution unless
                                      any conda option has been set explicitly.
      --galaxy_startup_timeout INTEGER RANGE
                                      Wait for galaxy to start before assuming
                                      Galaxy did not start.  [x>=1]
      --job_config_file FILE          Job configuration file for Galaxy to target.
      --tool_dependency_dir DIRECTORY
                                      Tool dependency dir for Galaxy to target.
      --tool_data_path DIRECTORY      Directory where data used by tools is located.
                                      Required if tests are run in docker and should
                                      make use of external reference data.
      --help                          Show this message and exit.
    
