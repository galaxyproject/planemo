
``rerun`` command
========================================

This section is auto-generated from the help text for the planemo command
``rerun``. This help message can be generated with ``planemo rerun
--help``.

**Usage**::

    planemo rerun [OPTIONS] RERUNNABLE_IDS

**Help**


Planemo command for rerunning failed jobs on an external Galaxy server.
Supply a list of history, invocation or job IDs, identifying the ID type using the
--invocation, --history or --job flag.

When passing the --invocation flag, the associated workflow will be rerun.
When using `--use_cache` (the default), successful jobs will be reused when available, so only failed jobs will actually
be executed on the compute node(s).

When not rerunning invocations, please note that attempting to rerun non-remappable jobs will result in an exit code of 1.
As jobs cannot be remapped more than once, running this command two or more times with the same
history or job IDs will therefore return an exit code of 1. If avoiding this is important,
you should specify the invocation ID instead if possible.

::

    % planemo rerun --invocation / --history / --job RERUNNABLE_IDS

**Options**::


      --profile TEXT                Name of profile (created with the profile_create
                                    command) to use with this command.
      --galaxy_url TEXT             Remote Galaxy URL to use with external Galaxy
                                    engine.
      --galaxy_user_key TEXT        User key to use with external Galaxy engine.
      --invocation                  Rerun failed jobs associated by one or more
                                    invocation IDs.
      --history                     Rerun failed jobs associated by one or more
                                    history IDs.
      --job                         Rerun failed jobs specified by one or more job
                                    IDs.
      --use_cache / --no_use_cache  Use cached job results if available.
      --help                        Show this message and exit.
    
