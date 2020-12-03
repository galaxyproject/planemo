
``workflow_test_init`` command
======================================

This section is auto-generated from the help text for the planemo command
``workflow_test_init``. This help message can be generated with ``planemo workflow_test_init
--help``.

**Usage**::

    planemo workflow_test_init [OPTIONS] WORKFLOW_PATH_OR_ID

**Help**

Initialize a Galaxy workflow test description for supplied workflow.

Be sure to your lint your workflow with ``workflow_lint`` before calling this
to ensure inputs and outputs comply with best practices that make workflow
testing easier.

**Options**::


      -f, --force                     Overwrite existing files if present.
      -o, --output FILE
      --split_test / --no_split_test  Write workflow job and test definitions to
                                      separate files.
    
      --help                          Show this message and exit.
    
