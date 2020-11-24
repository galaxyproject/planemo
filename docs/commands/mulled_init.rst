
``mulled_init`` command
======================================

This section is auto-generated from the help text for the planemo command
``mulled_init``. This help message can be generated with ``planemo mulled_init
--help``.

**Usage**::

    planemo mulled_init [OPTIONS]

**Help**

Download and install involucro for mull command.

This will happen automatically when using the mull command, but this can
be pre-installed in an environment using this command.

**Options**::


      --mulled_conda_version TEXT  Install a specific version of Conda before
                                   running the command, by default the version that
                                   comes with the continuumio miniconda3 image will
                                   be used under Linux and under Mac OS X Conda will
                                   be upgraded to to work around a bug in 4.2.
    
      --mulled_namespace TEXT      Build a mulled image with the specified namespace
                                   - defaults to biocontainers. Galaxy currently
                                   only recognizes images with the namespace
                                   biocontainers.
    
      --mulled_command TEXT        Mulled action to perform for targets - this
                                   defaults to 'build-and-test'.
    
      --help                       Show this message and exit.
    
