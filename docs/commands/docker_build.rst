
``docker_build`` command
======================================

This section is auto-generated from the help text for the planemo command
``docker_build``. This help message can be generated with ``planemo docker_build
--help``.

**Usage**::

    planemo docker_build [OPTIONS] TOOL_PATH

**Help**

Build (and optionally cache) Docker images.

Loads the tool or tools referenced by ``TOOL_PATH`` (by default all tools
in current directory), and ensures they all reference the same Docker image
and then attempts to find a Dockerfile for these tools (can be explicitly
specified with ``--dockerfile`` but by default it will check the tool's
directory and the current directory as well).

This command will then build and tag the image so it is ready to be tested
and published. The docker_shell command be used to test out the built
image.

::

    % planemo docker_build bowtie2.xml # asssumes Dockerfile in same dir
    % planemo docker_shell --from_tag bowtie2.xml

This can optionally also cache the images.

**Options**::


      --dockerfile TEXT
      --docker_image_cache TEXT
      --docker_cmd TEXT               Command used to launch docker (defaults to
                                      docker).
    
      --docker_sudo / --no_docker_sudo
                                      Flag to use sudo when running docker.
      --docker_sudo_cmd TEXT          sudo command to use when --docker_sudo is
                                      enabled (defaults to sudo).
    
      --docker_host TEXT              Docker host to target when executing docker
                                      commands (defaults to localhost).
    
      --help                          Show this message and exit.
    
