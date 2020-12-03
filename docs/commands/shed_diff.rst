
``shed_diff`` command
======================================

This section is auto-generated from the help text for the planemo command
``shed_diff``. This help message can be generated with ``planemo shed_diff
--help``.

**Usage**::

    planemo shed_diff [OPTIONS] PROJECT

**Help**

diff between local repository and Tool Shed.

By default, this will produce a diff between this repository and what
would be uploaded to the Tool Shed with the `shed_upload` command - but
this command can be made to compare other combinations of repositories.
Here are some examples

::

    $ # diff for this repository and the main Tool Shed
    $ planemo shed_diff
    $ # diff for this repository and the test Tool Shed
    $ planemo shed_diff --shed_target testtoolshed
    $ # diff for the test Tool Shed and main Tool Shed
    $ planemo shed_diff --shed_target_source testtoolshed
    $ # diff for two an explicitly specified repositories (ignores
    $ # current project's shed YAML file.)
    $ planemo shed_diff --owner peterjc --name blast_rbh
        --shed_target_source testtoolshed

This command will return an exit code of:

- 0 if there are no detected differences.
- 1 if there are differences.
- 2 if the target repository doesn't exist.
- >200 if there are errors attempting to perform a diff.

**Warning:** ``shed_diff`` doesn't inspect repository metadata, this
difference applies only to the file contents of files that would actually be
uploaded to the repository.

**Options**::


      -r, --recursive            Recursively perform command for nested repository
                                 directories.
    
      --fail_fast                If multiple repositories are specified and an error
                                 occurs stop immediately instead of processing
                                 remaining repositories.
    
      --owner TEXT               Tool Shed repository owner (username).
      --name TEXT                Tool Shed repository name (defaults to the inferred
                                 tool directory name).
    
      --shed_email TEXT          E-mail for Tool Shed auth (required unless shed_key
                                 is specified).
    
      --shed_key TEXT            API key for Tool Shed access. An API key is
                                 required unless e-mail and password is specified.
                                 This key can be specified with either --shed_key or
                                 --shed_key_from_env.
    
      --shed_key_from_env TEXT   Environment variable to read API key for Tool Shed
                                 access from.
    
      --shed_password TEXT       Password for Tool Shed auth (required unless
                                 shed_key is specified).
    
      -t, --shed_target TEXT     Tool Shed to target (this can be 'toolshed',
                                 'testtoolshed', 'local' (alias for
                                 http://localhost:9009/), an arbitrary url or
                                 mappings defined ~/.planemo.yml.
    
      -o, --output PATH          Send diff output to specified file.
      --shed_target_source TEXT  Source Tool Shed to diff against (will ignore local
                                 project info specified). To compare the main Tool
                                 Shed against the test, set this to testtoolshed.
    
      --raw                      Do not attempt smart diff of XML to filter out
                                 attributes populated by the Tool Shed.
    
      --report_xunit PATH        Output an XUnit report, useful for CI testing
      --help                     Show this message and exit.
    
