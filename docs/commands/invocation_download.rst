
``invocation_download`` command
========================================

This section is auto-generated from the help text for the planemo command
``invocation_download``. This help message can be generated with ``planemo invocation_download
--help``.

**Usage**::

    planemo invocation_download [OPTIONS] INVOCATION_ID

**Help**

Download output files from a completed Galaxy workflow invocation.

This command allows downloading output files after a workflow has been executed
through Galaxy's web interface or through planemo.

::

    % planemo invocation_download INVOCATION_ID

**Options**::


      --output_directory, --outdir DIRECTORY
                                      Where to store outputs of a 'run' task.
      --output_json FILE              Where to store JSON dictionary describing
                                      outputs of a 'run' task.
      --profile TEXT                  Name of profile (created with the
                                      profile_create command) to use with this
                                      command.
      --galaxy_url TEXT               Remote Galaxy URL to use with external Galaxy
                                      engine.
      --galaxy_user_key TEXT          User key to use with external Galaxy engine.
      --ignore_missing_output         Ignore missing output files
      --help                          Show this message and exit.
    
