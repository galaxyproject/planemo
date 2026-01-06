
``shed_lint`` command
========================================

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
      --shed_fail_fast                If multiple repositories are specified and an
                                      error occurs stop immediately instead of
                                      processing remaining repositories.
      --report_level [all|warn|error]
      --fail_level [warn|error]
      -s, --skip TEXT                 Comma-separated list of lint tests to skip
                                      (e.g. passing --skip 'citations,xml_order'
                                      would skip linting of citations and best-
                                      practice XML ordering.
      --skip_file FILE                File containing a list of lint tests to skip
      --tools                         Lint tools discovered in the process of
                                      linting repositories.
      --ensure_metadata               Ensure .shed.yml files contain enough metadata
                                      for each repository to allow automated
                                      creation and/or updates.
      --urls                          Check validity of URLs in XML files
      --biocontainer, --biocontainers
                                      Check best practice BioContainer namespaces
                                      for a container definition applicable for this
                                      tool.
      --help                          Show this message and exit.
    
