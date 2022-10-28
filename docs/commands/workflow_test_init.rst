
``workflow_test_init`` command
======================================

This section is auto-generated from the help text for the planemo command
``workflow_test_init``. This help message can be generated with ``planemo workflow_test_init
--help``.

**Usage**::

    planemo workflow_test_init [OPTIONS] WORKFLOW_PATH_OR_ID

**Help**

Initialize a Galaxy workflow test description for supplied workflow.

Be sure to lint your workflow with ``workflow_lint`` before calling this
to ensure inputs and outputs comply with best practices that make workflow
testing easier.

**Options**::


      -f, --force                     Overwrite existing files if present.
      -o, --output FILE
      --split_test / --no_split_test  Write workflow job and test definitions to
                                      separate files.
      --galaxy_url TEXT               Remote Galaxy URL to use with external Galaxy
                                      engine.
      --galaxy_user_key TEXT          User key to use with external Galaxy engine.
      --from_invocation / --from_uri  Build a workflow test or job description from
                                      an invocation ID run on an external Galaxy.A
                                      Galaxy URL and API key must also be specified.
                                      This allows test data to be downloadedand
                                      inputs and parameters defined automatically.
                                      Alternatively, the default is to build
                                      thedescriptions from a provided workflow URI.
      --profile TEXT                  Name of profile (created with the
                                      profile_create command) to use with this
                                      command.
      --help                          Show this message and exit.
    
