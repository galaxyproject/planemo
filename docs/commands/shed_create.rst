
``shed_create`` command
======================================

This section is auto-generated from the help text for the planemo command
``shed_create``. This help message can be generated with ``planemo shed_create
--help``.

**Usage**::

    planemo shed_create [OPTIONS] PROJECT

**Help**

Create a repository in a Galaxy Tool Shed.

This will read the settings from the ``.shed.yml`` file.

**Options**::


      -r, --recursive           Recursively perform command for nested repository
                                directories.
    
      --fail_fast               If multiple repositories are specified and an error
                                occurs stop immediately instead of processing
                                remaining repositories.
    
      --owner TEXT              Tool Shed repository owner (username).
      --name TEXT               Tool Shed repository name (defaults to the inferred
                                tool directory name).
    
      --shed_email TEXT         E-mail for Tool Shed auth (required unless shed_key
                                is specified).
    
      --shed_key TEXT           API key for Tool Shed access. An API key is required
                                unless e-mail and password is specified. This key
                                can be specified with either --shed_key or
                                --shed_key_from_env.
    
      --shed_key_from_env TEXT  Environment variable to read API key for Tool Shed
                                access from.
    
      --shed_password TEXT      Password for Tool Shed auth (required unless
                                shed_key is specified).
    
      -t, --shed_target TEXT    Tool Shed to target (this can be 'toolshed',
                                'testtoolshed', 'local' (alias for
                                http://localhost:9009/), an arbitrary url or
                                mappings defined ~/.planemo.yml.
    
      -m, --message TEXT        Commit message for tool shed upload.
      --skip_upload             Skip upload contents as part of operation, only
                                update metadata.
    
      --help                    Show this message and exit.
    
