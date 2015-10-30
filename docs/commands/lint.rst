
``lint`` command
======================================

This section is auto-generated from the help text for the planemo command
``lint``. This help message can be generated with ``planemo lint
--help``.

**Usage**::

    planemo lint [OPTIONS] TOOL_PATH

**Help**

Check specified tool(s) for common errors and adherence to best
practices.

**Options**::


      --report_level [all|warn|error]
      --report_xunit PATH             Output an XUnit report, useful for CI
                                      testing
      --fail_level [warn|error]
      -s, --skip TEXT                 Comma-separated list of lint tests to skip
                                      (e.g send .--skip 'citations,xml_order' to
                                      skip linting of citations and best-practice
                                      XML ordering.
      --xsd                           Include experimental tool XSD validation in
                                      linting process (requires xmllint on PATH or
                                      lxml installed).
      -r, --recursive                 Recursively perform command for
                                      subdirectories.
      --urls                          Check validity of URLs in XML files
      --help                          Show this message and exit.
    
