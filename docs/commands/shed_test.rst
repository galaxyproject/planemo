
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
      --fail_fast                     If multiple repositories are specified and
                                      an error occurs stop immediately instead of
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
      --shed_key_from_env TEXT        Environment variable to read API key for
                                      Tool Shed access from.
      --shed_password TEXT            Password for Tool Shed auth (required unless
                                      shed_key is specified).
      -t, --shed_target TEXT          Tool Shed to target (this can be 'toolshed',
                                      'testtoolshed', 'local' (alias for
                                      http://localhost:9009/), an arbitrary url or
                                      mappings defined ~/.planemo.yml.
      --galaxy_root DIRECTORY         Root of development galaxy directory to
                                      execute command with.
      --install_galaxy                Download and configure a disposable copy of
                                      Galaxy from github.
      --no_cache_galaxy               Skip caching of Galaxy source and
                                      dependencies obtained with --install_galaxy.
                                      Not caching this results in faster downloads
                                      (no git) - so is better on throw away
                                      instances such with TravisCI.
      --no_cleanup                    Do not cleanup temp files created for and by
                                      Galaxy.
      --job_config_file PATH          Job configuration file for Galaxy to target.
      --update_test_data              Update test-data directory with job outputs
                                      (normally written to directory
                                      --job_output_files if specified.)
      --test_output PATH              Output test report (HTML - for humans)
                                      defaults to tool_test_output.html.
      --test_output_text PATH         Output test report (Basic text - for display
                                      in CI)
      --test_output_markdown PATH     Output test report (Markdown style - for
                                      humans & computers)
      --test_output_xunit PATH        Output test report (xUnit style - for
                                      computers).
      --test_output_json PATH         Output test report (planemo json) defaults
                                      to tool_test_output.json.
      --job_output_files DIRECTORY    Write job outputs to specified directory.
      --summary [none|minimal|compact]
                                      Summary style printed to planemo's standard
                                      output (see output reports for more complete
                                      summary). Set to 'none' to disable
                                      completely.
      --skip_dependencies             Do not install shed dependencies as part of
                                      repository installation.
      --help                          Show this message and exit.
    
