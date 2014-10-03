
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


      --test_output PATH       Output test report.
      --galaxy_root DIRECTORY  Root of development galaxy directory to execute
                               command with.
      --install_galaxy         Download and configure a disposable copy of Galaxy
                               from github.
      --test-data DIRECTORY    test-data directory to for specified tool(s).
      --help                   Show this message and exit.
    
