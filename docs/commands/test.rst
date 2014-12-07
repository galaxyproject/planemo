
``test`` command
===============================

This section is auto-generated from the help text for the planemo command
``test``. This help message can be generated with ``planemo test
--help``.

**Usage**::

    planemo test [OPTIONS] TOOL_PATH

**Help**

Run the tests in the specified tool tests in a Galaxy instance.

All referenced tools (by default all the tools in the current working
directory) will be tested and the resulted disposited in path specified
with ``--test_output`` (defaults to tool_test_output.html).

To run these tests planemo needs a Galaxy instance to utilize, planemo
will search parent directories to see if any is a Galaxy instance
- but one can pick the Galaxy instance to use with the --galaxy_root
option or force planemo to download a disposable instance with the
``--install_galaxy`` flag.

planemo uses temporarily generated config files and environment variables
to attempt to shield this execution of Galaxy from manually launched runs
against that same Galaxy root - but this may not be bullet proof yet so
please careful and do not try this against production Galaxy instances.

**Options**::


      --test_output PATH              Output test report (HTML - for humans).
      --test_output_xunit PATH        Output test report (xUnit style - for
                                      computers).
      --job_output_files DIRECTORY    Write job outputs to specified directory.
      --update_test_data              Update test-data directory with job outputs
                                      (normally written to directory
                                      --job_output_files if specified.)
      --summary [none|minimal|compact]
                                      Summary style printed to planemo's standard
                                      output (see output reports for more complete
                                      summary). Set to 'none' to disable
                                      completely.
      --galaxy_root DIRECTORY         Root of development galaxy directory to
                                      execute command with.
      --install_galaxy                Download and configure a disposable copy of
                                      Galaxy from github.
      --test-data DIRECTORY           test-data directory to for specified
                                      tool(s).
      --dependency_resolvers_config_file PATH
                                      Dependency resolver configuration for Galaxy
                                      to target.
      --job_config_file PATH          Job configuration file for Galaxy to target.
      --tool_dependency_dir DIRECTORY
                                      Tool dependency dir for Galaxy to target.
      --brew_dependency_resolution    Configure Galaxy to use brew dependency
                                      resolution.
      --help                          Show this message and exit.
    
