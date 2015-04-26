
``shed_create`` command
======================================

This section is auto-generated from the help text for the planemo command
``shed_create``. This help message can be generated with ``planemo shed_create
--help``.

**Usage**::

    planemo shed_create [OPTIONS] PROJECT

**Help**

Create a repository in a Galaxy Tool Shed from a ``.shed.yml`` file.

**Options**::


      --owner TEXT          Tool Shed repository owner (username).
      --name TEXT           Tool Shed repository name (defaults to the inferred
                            tool directory name).
      --shed_target TEXT    Tool Shed to target (this can be 'toolshed',
                            'testtoolshed', 'local' (alias for
                            http://localhost:9009/) or an arbitraryurl).
      --shed_key TEXT       API key for Tool Shed access (required unless
                            e-mail/pass specified).
      --shed_email TEXT     E-mail for Tool Shed auth (required unless shed_key is
                            specified).
      --shed_password TEXT  Password for Tool Shed auth (required unless shed_key
                            is specified).
      -r, --recursive       Recursively perform command for nested repository
                            directories.
      --help                Show this message and exit.
    
