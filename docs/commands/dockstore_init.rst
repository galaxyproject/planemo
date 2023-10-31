
``dockstore_init`` command
========================================

This section is auto-generated from the help text for the planemo command
``dockstore_init``. This help message can be generated with ``planemo dockstore_init
--help``.

**Usage**::

    planemo dockstore_init [OPTIONS] PROJECT

**Help**

Initialize a .dockstore.yml configuration file for workflows in directory.

Walk supplied directory and find all Galaxy workflows and test configurations
and create a ``.dockstore.yml`` with references to these files. Be sure to push
this file to Github before registering your workflow repository with Dockstore.

Visit Dockstore at https://dockstore.org/. Find more about registering workflows
with Dockstore at
https://docs.dockstore.org/en/develop/getting-started/dockstore-workflows.html.

**Options**::


      --publish / --no_publish  Set publish attribute to true in .dockstore.yml file
      --help                    Show this message and exit.
    
