
``shed_lint`` command
======================================

This section is auto-generated from the help text for the planemo command
``shed_lint``. This help message can be generated with ``planemo shed_lint
--help``.

**Usage**::

    planemo shed_lint [OPTIONS] PROJECT

**Help**

Check a Tool Shed repository for common problems.

**Options**::


      -r, --recursive                 Recursively perform command for nested
                                      repository directories.
      --fail_fast                     If multiple repositories are specified and
                                      an error occurs stop immediately instead of
                                      processing remaining repositories.
      --report_level [all|warn|error]
      --fail_level [warn|error]
      --tools                         Lint tools discovered in the process of
                                      linting repositories.
      --xsd                           Include experimental tool XSD validation in
                                      linting process (requires xmllint on PATH or
                                      lxml installed).
      --ensure_metadata               Ensure .shed.yml files contain enough
                                      metadata for each repository to allow
                                      automated creation and/or updates.
      --help                          Show this message and exit.
    
