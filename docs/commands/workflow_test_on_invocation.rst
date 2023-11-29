
``workflow_test_on_invocation`` command
========================================

This section is auto-generated from the help text for the planemo command
``workflow_test_on_invocation``. This help message can be generated with ``planemo workflow_test_on_invocation
--help``.

**Usage**::

    planemo workflow_test_on_invocation [OPTIONS] TEST DEFINITION INVOCATION

**Help**

Run defined tests against existing workflow invocation.
**Options**::


      --galaxy_url TEXT               Remote Galaxy URL to use with external Galaxy
                                      engine.  [required]
      --galaxy_user_key TEXT          User key to use with external Galaxy engine.
                                      [required]
      --test_index INTEGER            Select which test to check. Counting starts at
                                      1
      --update_test_data              Update test-data directory with job outputs
                                      (normally written to directory
                                      --job_output_files if specified.)
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
      --test_output_allure DIRECTORY  Output test allure2 framework resutls
      --test_output_json PATH         Output test report (planemo json) defaults to
                                      tool_test_output.json.
      --job_output_files DIRECTORY    Write job outputs to specified directory.
      --summary [none|minimal|compact]
                                      Summary style printed to planemo's standard
                                      output (see output reports for more complete
                                      summary). Set to 'none' to disable completely.
      --test_timeout INTEGER          Maximum runtime of a single test in seconds.
      --help                          Show this message and exit.
    
