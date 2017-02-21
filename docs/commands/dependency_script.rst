
``dependency_script`` command
======================================

This section is auto-generated from the help text for the planemo command
``dependency_script``. This help message can be generated with ``planemo dependency_script
--help``.

**Usage**::

    planemo dependency_script [OPTIONS] PROJECT

**Help**

Compile tool_dependencies.xml to bash script.

An experimental approach parsing tool_dependencies.xml files into
bash shell scripts, intended initially for use within Continuous
Integration testing setups like TravisCI.

Parses the ``tool_dependencies.xml`` files from the specified projects,
and converts them into an installation bash script (``dep_install.sh``),
and a shell script (``env.sh``) defining any new/edited environment
variables.

These are intended to be used via ``bash dep_install.sh`` (once), and as
``source env.sh`` prior to running any of the dependencies to set the
environment variable within the current shell session.

Both ``dep_install.sh`` and ``env.sh`` require ``$INSTALL_DIR`` be defined
before running them, set to an existing directory with write permissions.
Beware than if run on multiple tools, they can over-write each other (for
example if you have packages for different versions of the same tool). In
this case make separate calls to ``planemo dependency_script`` and call
the scripts with different installation directories.

This command will download (and cache) any URLs specified via Galaxy
download actions. This is in order to decompress them and determine the
relevant sub-folder to change into as per the Tool Shed install mechanism,
so that this can be recorded as a ``cd`` comand in the bash script.

The download cache used by ``planemo dependency_script`` and the resulting
output script ``dep_install.sh`` defaults to ``./download_cache`` (under
the current working directory), and can be set with ``$DOWNLOAD_CACHE``.

If the ``tool_dependencies.xml`` file includes SHA256 checksums for
downloads, these will be verified after downloading to the cache (by
either ``planemo dependency_script`` or ``bash dep_install.sh``).

This is experimental, and is initially intended for use within continuous
integration testing setups like TravisCI to both verify the dependency
installation receipe works, and to use this to run functional tests.

**Options**::


      -r, --recursive             Recursively perform command for nested repository
                                  directories.
      --fail_fast                 If multiple repositories are specified and an
                                  error occurs stop immediately instead of
                                  processing remaining repositories.
      --download_cache DIRECTORY  Directory to cache downloaded files, default is
                                  $DOWNLOAD_CACHE
      --help                      Show this message and exit.
    
