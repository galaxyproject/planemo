
``lint`` command
======================================

This section is auto-generated from the help text for the planemo command
``lint``. This help message can be generated with ``planemo lint
--help``.

**Usage**::

    planemo lint [OPTIONS] TOOL_PATH

**Help**

Check for common errors and best practices.
**Options**::


      --report_level [all|warn|error]
      --report_xunit PATH             Output an XUnit report, useful for CI testing
      --fail_level [warn|error]
      -s, --skip TEXT                 Comma-separated list of lint tests to skip
                                      (e.g. passing --skip 'citations,xml_order'
                                      would skip linting of citations and best-
                                      practice XML ordering.
    
      --xsd / --no_xsd                Include tool XSD validation in linting
                                      process.
    
      -r, --recursive                 Recursively perform command for
                                      subdirectories.
    
      --urls                          Check validity of URLs in XML files
      --doi                           Check validity of DOIs in XML files
      --conda_requirements            Check tool requirements for availability in
                                      best practice Conda channels.
    
      --biocontainer, --biocontainers
                                      Check best practice BioContainer namespaces
                                      for a container definition applicable for this
                                      tool.
    
      --help                          Show this message and exit.
    
