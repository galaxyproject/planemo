
``test`` command
======================================

This section is auto-generated from the help text for the planemo command
``test``. This help message can be generated with ``planemo test
--help``.

**Usage**::

    planemo test [OPTIONS] TOOL_PATH

**Help**

Run the tests in the specified tool tests in a Galaxy instance.

All referenced tools (by default all the tools in the current working
directory) will be tested and the results quickly summarized.

To run these tests planemo needs a Galaxy instance to utilize, planemo
will search parent directories to see if any is a Galaxy instance
- but one can pick the Galaxy instance to use with the --galaxy_root
option or force planemo to download a disposable instance with the
``--install_galaxy`` flag.

In additon to to quick summary printed to the console - various detailed
output summaries can be configured. ``tool_test_output.html`` (settable
via ``--test_output``) will contain a human consumable HTML report
describing the test run. A JSON file (settable via ``--test_output_json``
and defaulting to ``tool_test_output.json``) will also be created. These
files can can be disabled by passing in empty arguments or globally by
setting the values ``default_test_output`` and/or
``default_test_output_json`` in ``~/.planemo.yml`` to ``null``. For
continuous integration testing a xUnit-style report can be confiured using
the ``--test_output_xunit``.

planemo uses temporarily generated config files and environment variables
to attempt to shield this execution of Galaxy from manually launched runs
against that same Galaxy root - but this may not be bullet proof yet so
please careful and do not try this against production Galaxy instances.

**Options**::


      --failed                        Re-run only failed tests. This command will
                                      read tool_test_output.json to determine
                                      which tests failed so this file must have
                                      been produced with the same set of tool ids
                                      previously.
      --galaxy_root DIRECTORY         Root of development galaxy directory to
                                      execute command with.
      --galaxy_sqlite_database DIRECTORY
                                      Preseeded Galaxy sqlite database to target.
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
      --test_data DIRECTORY           test-data directory to for specified
                                      tool(s).
      --tool_data_table PATH          tool_data_table_conf.xml file to for
                                      specified tool(s).
      --dependency_resolvers_config_file PATH
                                      Dependency resolver configuration for Galaxy
                                      to target.
      --tool_dependency_dir DIRECTORY
                                      Tool dependency dir for Galaxy to target.
      --brew_dependency_resolution    Configure Galaxy to use plain brew
                                      dependency resolution.
      --shed_dependency_resolution    Configure Galaxy to use brewed Tool Shed
                                      dependency resolution.
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
      --help                          Show this message and exit.
    
