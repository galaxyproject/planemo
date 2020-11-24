
``shed_test`` command
======================================

This section is auto-generated from the help text for the planemo command
``shed_test``. This help message can be generated with ``planemo shed_test
--help``.

**Usage**::

    planemo shed_test [OPTIONS] PROJECT

**Help**

Run tests of published shed artifacts.

This command will start a Galaxy instance configured to target the
specified shed, find published artifacts (tools and dependencies)
corresponding to command-line arguments and ``.shed.yml`` file(s),
install these artifacts, and run the tool tests for these commands.

This command requires the target to be version 15.07 or newer.

**Options**::


      -r, --recursive                 Recursively perform command for nested
                                      repository directories.
    
      --fail_fast                     If multiple repositories are specified and an
                                      error occurs stop immediately instead of
                                      processing remaining repositories.
    
      --owner TEXT                    Tool Shed repository owner (username).
      --name TEXT                     Tool Shed repository name (defaults to the
                                      inferred tool directory name).
    
      --shed_email TEXT               E-mail for Tool Shed auth (required unless
                                      shed_key is specified).
    
      --shed_key TEXT                 API key for Tool Shed access. An API key is
                                      required unless e-mail and password is
                                      specified. This key can be specified with
                                      either --shed_key or --shed_key_from_env.
    
      --shed_key_from_env TEXT        Environment variable to read API key for Tool
                                      Shed access from.
    
      --shed_password TEXT            Password for Tool Shed auth (required unless
                                      shed_key is specified).
    
      -t, --shed_target TEXT          Tool Shed to target (this can be 'toolshed',
                                      'testtoolshed', 'local' (alias for
                                      http://localhost:9009/), an arbitrary url or
                                      mappings defined ~/.planemo.yml.
    
      --galaxy_root DIRECTORY         Root of development galaxy directory to
                                      execute command with.
    
      --galaxy_python_version [3|3.6|3.7|3.8|3.9]
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
                                      --docker).
    
      --job_config_file FILE          Job configuration file for Galaxy to target.
      --tool_dependency_dir DIRECTORY
                                      Tool dependency dir for Galaxy to target.
      --update_test_data              Update test-data directory with job outputs
                                      (normally written to directory
                                      --job_output_files if specified.)
    
      --paste_test_data_paths / --no_paste_test_data_paths
                                      By default Planemo will use or not use
                                      Galaxy's path paste option to load test data
                                      into a history based on the engine type it is
                                      targeting. This can override the logic to
                                      explicitly enable or disable path pasting.
    
      --test_output PATH              Output test report (HTML - for humans)
                                      defaults to tool_test_output.html.
    
      --test_output_text PATH         Output test report (Basic text - for display
                                      in CI)
    
      --test_output_markdown PATH     Output test report (Markdown style - for
                                      humans & computers)
    
      --test_output_xunit PATH        Output test report (xunit style - for CI
                                      systems
    
      --test_output_junit PATH        Output test report (jUnit style - for CI
                                      systems
    
      --test_output_json PATH         Output test report (planemo json) defaults to
                                      tool_test_output.json.
    
      --job_output_files DIRECTORY    Write job outputs to specified directory.
      --summary [none|minimal|compact]
                                      Summary style printed to planemo's standard
                                      output (see output reports for more complete
                                      summary). Set to 'none' to disable completely.
    
      --skip_dependencies             Do not install shed dependencies as part of
                                      repository installation.
    
      --help                          Show this message and exit.
    
