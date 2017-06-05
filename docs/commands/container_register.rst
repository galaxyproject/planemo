
``container_register`` command
======================================

This section is auto-generated from the help text for the planemo command
``container_register``. This help message can be generated with ``planemo container_register
--help``.

**Usage**::

    planemo container_register [OPTIONS] TOOL_PATH

**Help**

Register multi-requirement containers as needed.

BioContainers publishes all Bioconda packages automatically as individual
container images. These however are not enough for tools with multiple
best-practice requirements. Such requirements should be recorded and published
so that a container can be created and registered for these tools.

**Options**::


      -r, --recursive                 Recursively perform command for
                                      subdirectories.
      --mulled_namespace TEXT         Build a mulled image with the specified
                                      namespace - defaults to biocontainers. Galaxy
                                      currently only recognizes images with the
                                      namespace biocontainers.
      --output_directory DIRECTORY    Container registration directory (defaults to
                                      ~/.planemo/multireqcontainers.
      -m, --message TEXT              Commit and pull request message template for
                                      registration interactions.
      --pull_request / --no_pull_request
                                      Fork and create a pull request against
                                      jmchilton/multireqcontainers for these
                                      changes.
      --force_push / --no_force_push  Force push branch for pull request in case it
                                      already exists.
      --help                          Show this message and exit.
    
