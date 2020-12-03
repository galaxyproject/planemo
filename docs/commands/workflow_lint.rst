
``workflow_lint`` command
======================================

This section is auto-generated from the help text for the planemo command
``workflow_lint``. This help message can be generated with ``planemo workflow_lint
--help``.

**Usage**::

    planemo workflow_lint [OPTIONS] TARGET

**Help**

Check workflows for syntax errors and best practices.
**Options**::


      --report_level [all|warn|error]
      --report_xunit PATH             Output an XUnit report, useful for CI testing
      --fail_level [warn|error]
      -s, --skip TEXT                 Comma-separated list of lint tests to skip
                                      (e.g. passing --skip 'citations,xml_order'
                                      would skip linting of citations and best-
                                      practice XML ordering.
    
      --help                          Show this message and exit.
    
