
``shed_upload`` command
======================================

This section is auto-generated from the help text for the planemo command
``shed_upload``. This help message can be generated with ``planemo shed_upload
--help``.

**Usage**::

    planemo shed_upload [OPTIONS] PROJECT

**Help**

Low-level command to upload tarballs.

Generally, ``shed_update`` should be used instead since it also updates
both tool shed contents (via tar ball generation and upload) as well as
metadata (to handle metadata changes in ``.shed.yml`` files).

::

    % planemo shed_upload --tar_only  ~/
    % tar -tzf shed_upload.tar.gz
    test-data/blastdb.loc
    ...
    tools/ncbi_blast_plus/tool_dependencies.xml
    % tar -tzf shed_upload.tar.gz | wc -l
    117


**Options**::


      -r, --recursive              Recursively perform command for nested repository
                                   directories.
    
      --fail_fast                  If multiple repositories are specified and an
                                   error occurs stop immediately instead of
                                   processing remaining repositories.
    
      --owner TEXT                 Tool Shed repository owner (username).
      --name TEXT                  Tool Shed repository name (defaults to the
                                   inferred tool directory name).
    
      --shed_email TEXT            E-mail for Tool Shed auth (required unless
                                   shed_key is specified).
    
      --shed_key TEXT              API key for Tool Shed access. An API key is
                                   required unless e-mail and password is specified.
                                   This key can be specified with either --shed_key
                                   or --shed_key_from_env.
    
      --shed_key_from_env TEXT     Environment variable to read API key for Tool
                                   Shed access from.
    
      --shed_password TEXT         Password for Tool Shed auth (required unless
                                   shed_key is specified).
    
      -t, --shed_target TEXT       Tool Shed to target (this can be 'toolshed',
                                   'testtoolshed', 'local' (alias for
                                   http://localhost:9009/), an arbitrary url or
                                   mappings defined ~/.planemo.yml.
    
      -m, --message TEXT           Commit message for tool shed upload.
      --force_repository_creation  If a repository cannot be found for the specified
                                   user/repo name pair, then automatically create
                                   the repository in the toolshed.
    
      --check_diff                 Skip uploading if the shed_diff detects there
                                   would be no 'difference' (only attributes
                                   populated by the shed would be updated.)
    
      --tar_only                   Produce tar file for upload but do not publish to
                                   a tool shed.
    
      --tar FILE                   Specify a pre-existing tar file instead of
                                   automatically building one as part of this
                                   command.
    
      --help                       Show this message and exit.
    
