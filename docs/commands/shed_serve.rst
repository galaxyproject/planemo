
``shed_serve`` command
======================================

This section is auto-generated from the help text for the planemo command
``shed_serve``. This help message can be generated with ``planemo shed_serve
--help``.

**Usage**::

    planemo shed_serve [OPTIONS] PROJECT

**Help**

 Serve a transient Galaxy instance after installing repositories
from a remote Tool Shed.

**Options**::


      -r, --recursive          Recursively perform command for nested repository
                               directories.
      --fail_fast              If multiple repositories are specified and an error
                               occurs stop immediately instead of processing
                               remaining repositories.
      --owner TEXT             Tool Shed repository owner (username).
      --name TEXT              Tool Shed repository name (defaults to the inferred
                               tool directory name).
      --shed_email TEXT        E-mail for Tool Shed auth (required unless shed_key
                               is specified).
      --shed_key TEXT          API key for Tool Shed access (required unless
                               e-mail/pass specified).
      --shed_password TEXT     Password for Tool Shed auth (required unless
                               shed_key is specified).
      --shed_target TEXT       Tool Shed to target (this can be 'toolshed',
                               'testtoolshed', 'local' (alias for
                               http://localhost:9009/) or an arbitraryurl).
      --galaxy_root DIRECTORY  Root of development galaxy directory to execute
                               command with.
      --install_galaxy         Download and configure a disposable copy of Galaxy
                               from github.
      --no_cache_galaxy        Skip caching of Galaxy source and dependencies
                               obtained with --install_galaxy. Not caching this
                               results in faster downloads (no git) - so is better
                               on throw away instances such with TravisCI.
      --no_cleanup             Do not cleanup temp files created for and by
                               Galaxy.
      --job_config_file PATH   Job configuration file for Galaxy to target.
      --port INTEGER           Port to serve Galaxy on (default is 9090).
      --skip_dependencies      Do not install shed dependencies as part of
                               repository installation.
      --help                   Show this message and exit.
    
