
``shed_update`` command
======================================

This section is auto-generated from the help text for the planemo command
``shed_update``. This help message can be generated with ``planemo shed_update
--help``.

**Usage**::

    planemo shed_update [OPTIONS] PROJECT

**Help**

Update Tool Shed repository.

By default this command will update both repository metadata
from ``.shed.yml`` and upload new contents from the repository
directory.

::

    % planemo shed_update

This will update the main tool shed with the repository defined
by a ``.shed.yml`` file in the current working directory. Both
the location of the ``.shed.yml`` and the tool shed to upload to
can be easily configured. For instance, the following command can
be used if ``.shed.yml`` if contained in ``path/to/repo`` and the
desire is to update the test tool shed.

::

    % planemo shed_update --shed_target testtoolshed path/to/repo

Another important option is ``--check_diff`` - this doesn't affect the
updating of shed metadata but it will check for differences before
uploading new contents to the tool shed. This may important because the
tool shed will automatically populate certain attributes in tool shed
artifact files (such as ``tool_dependencies.xml``) and this may
cause unwanted installable revisions to be created when there are no
important changes.

The lower-level ``shed_upload`` command should be used instead if
the repository doesn't define complete metadata in a ``.shed.yml``.

**Options**::


      --report_xunit PATH          Output an XUnit report, useful for CI testing
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
    
      --skip_upload                Skip upload contents as part of operation, only
                                   update metadata.
    
      --skip_metadata              Skip metadata update as part of operation, only
                                   upload new contents.
    
      --help                       Show this message and exit.
    
