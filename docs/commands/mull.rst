
``mull`` command
======================================

This section is auto-generated from the help text for the planemo command
``mull``. This help message can be generated with ``planemo mull
--help``.

**Usage**::

    planemo mull [OPTIONS] TOOL_PATH

**Help**

Build containers for specified tools.

Supplied tools will be inspected for referenced requirement packages. For
each combination of requirements a "mulled" container will be built. Galaxy
can automatically discover this container and subsequently use it to run
or test the tool.

For this to work, the tool's requirements will need to be present in a known
Conda channel such as bioconda (https://github.com/bioconda/bioconda-recipes).
This can be verified by running ``planemo lint --conda_requirements`` on the
target tool(s).

**Options**::


      -r, --recursive  Recursively perform command for subdirectories.
      --help           Show this message and exit.
    
