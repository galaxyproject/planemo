
``shed_upload`` command
===============================

This section is auto-generated from the help text for the planemo command
``shed_upload``. This help message can be generated with ``planemo shed_upload
--help``.

**Usage**::

    planemo shed_upload [OPTIONS] PROJECT

**Help**

Upload a tool directory as a tarball to a tool shed.

**Options**::


      --message TEXT        Commit message for tool shed upload.
      --owner TEXT          Tool shed repository owner (username).
      --name TEXT           Repository name (default to tool directory name).
      --shed_target TEXT    Tool shed to target (toolshed/testtoolshed/local/url).
      --shed_key TEXT       API key for shed access (required unless e-mail/pass
                            specified).
      --shed_email TEXT     E-mail for shed auth (required unless shed_key is
                            specified).
      --shed_password TEXT  Password for shed auth (required unless shed_key is
                            specified).
      --help                Show this message and exit.
    
