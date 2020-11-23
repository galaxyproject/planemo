
``shed_lint`` command
======================================

This section is auto-generated from the help text for the planemo command
``shed_lint``. This help message can be generated with ``planemo shed_lint
--help``.

**Usage**::

    planemo shed_lint [OPTIONS] PROJECT

**Help**

Check Tool Shed repository for common issues.

With the ``--tools`` flag, this command lints actual Galaxy tools
in addition to tool shed artifacts.

With the ``--urls`` flag, this command searches for
``<package>$URL</package>`` and download actions which specify URLs. Each
of those are accessed individually. By default, this tool requests the
first hundred or so bytes of each listed URL and validates that a 200 OK
was received. In tool XML files, the ``--urls`` option checks through the
help text for mentioned URLs and checks those.

**Options**::


      -r, --recursive                 Recursively perform command for nested
                                      repository directories.
    
      --fail_fast                     If multiple repositories are specified and an
                                      error occurs stop immediately instead of
                                      processing remaining repositories.
    
      --report_level [all|warn|error]
      --fail_level [warn|error]
      --tools                         Lint tools discovered in the process of
                                      linting repositories.
    
      --xsd / --no_xsd                Include tool XSD validation in linting
                                      process.
    
      --ensure_metadata               Ensure .shed.yml files contain enough metadata
                                      for each repository to allow automated
                                      creation and/or updates.
    
      --urls                          Check validity of URLs in XML files
      --help                          Show this message and exit.
    
