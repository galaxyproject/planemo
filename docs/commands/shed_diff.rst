
``shed_diff`` command
===============================

This section is auto-generated from the help text for the planemo command
``shed_diff``. This help message can be generated with ``planemo shed_diff
--help``.

**Usage**::

    planemo shed_diff [OPTIONS] PROJECT

**Help**

Produce diff between local repository and Tool Shed contents.

By default, this will produce a diff between this repository and what
would be uploaded to the Tool Shed with the `shed_upload` command - but
this command can be made to compare other combinations of repositories.
Here are some examples::

    % # diff for this repository and the main Tool Shed
    % planemo shed_diff
    % # diff for this repository and the test Tool Shed
    % planemo shed_diff --shed_target testtoolshed
    % # diff for the test Tool Shed and main Tool Shed
    % planemo shed_diff --shed_target_source testtoolshed
    % # diff for two an explicitly specified repositories (ignores
    % # current project's shed YAML file.)
    % planemo shed_diff --owner peterjc --name blast_rbh
        --shed_target_source testtoolshed

**Options**::


      --owner TEXT               Tool Shed repository owner (username).
      --name TEXT                Tool Shed repository name (defaults to the
                                 inferred tool directory name).
      --shed_target TEXT         Tool Shed to target (this can be 'toolshed',
                                 'testtoolshed', 'local' (alias for
                                 http://localhost:9009/) or an arbitraryurl).
      --shed_target_source TEXT  Source Tool Shed to diff against (will ignore
                                 local project info specified). To compare the
                                 main Tool Shed against the test, set this to
                                 testtoolshed.
      --help                     Show this message and exit.
    
