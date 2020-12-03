
``docker_shell`` command
======================================

This section is auto-generated from the help text for the planemo command
``docker_shell``. This help message can be generated with ``planemo docker_shell
--help``.

**Usage**::

    planemo docker_shell [OPTIONS] TOOL_PATH

**Help**

Launch shell in Docker container for a tool.

Will launch a shell in the Docker container referenced by the specified
tool. Prints a command to do this the way Galaxy would in job files it
generates - so be sure to wrap this in $(...) to launch the subshell.

::

    $ $(planemo docker_shell bowtie2.xml)
    ...
    root@b8754062f875:/#


**Options**::


      --from_tag                      Treat the tool's Docker container identifier
                                      as a locally cached tag.
    
      --shell TEXT                    Shell to launch in container (defaults to
                                      /bin/bash).
    
      --docker_cmd TEXT               Command used to launch docker (defaults to
                                      docker).
    
      --docker_sudo / --no_docker_sudo
                                      Flag to use sudo when running docker.
      --docker_sudo_cmd TEXT          sudo command to use when --docker_sudo is
                                      enabled (defaults to sudo).
    
      --docker_host TEXT              Docker host to target when executing docker
                                      commands (defaults to localhost).
    
      --help                          Show this message and exit.
    
