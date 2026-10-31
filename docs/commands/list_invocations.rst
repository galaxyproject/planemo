
``list_invocations`` command
========================================

This section is auto-generated from the help text for the planemo command
``list_invocations``. This help message can be generated with ``planemo list_invocations
--help``.

**Usage**::

    planemo list_invocations [OPTIONS] [WORKFLOW_IDENTIFIER]

**Help**


Get invocations, optionally filtering by workflow ID or alias.

**Options**::


      --raw                           Output invocations as JSON.
      --max-items, --max_items INTEGER RANGE
                                      Maximum number of invocations to return.
                                      [default: 100; x>=0]
      --offset-items, --offset_items INTEGER RANGE
                                      Skip this many invocations before returning
                                      results.  [default: 0; x>=0]
      --profile TEXT                  Name of profile (created with the
                                      profile_create command) to use with this
                                      command.  [required]
      --help                          Show this message and exit.
    
