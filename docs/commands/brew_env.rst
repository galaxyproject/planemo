
``brew_env`` command
===============================

This section is auto-generated from the help text for the planemo command
``brew_env``. This help message can be generated with ``planemo brew_env
--help``.

**Usage**::

    planemo brew_env [OPTIONS] TOOL_PATH

**Help**

Display commands used to modify environment to inject tool's brew
dependencies.::

    % . <(planemo brew_env bowtie2.xml)
    % which bowtie2
    /home/john/.linuxbrew/Cellar/bowtie2/2.1.0/bin/bowtie2

By default this will attempt to attempt to install these recipes as needed.
This automatic installation can be skipped with the ``--skip_install`` flag.

Intead of injecting the enviornment into your current shell using the above
idiom, the ``--shell`` flag can be sent to launch a new subshell when
sourced.::

    % . <(planemo brew_env --skip_install --shell bowtie2.xml)
    (bowtie2) % which bowtie2
    /home/john/.linuxbrew/Cellar/bowtie2/2.1.0/bin/bowtie2


**Options**::


      --brew PATH     Homebrew 'brew' executable to use.
      --skip_install  Skip installation - only source requirements already
                      available.
      --shell
      --help          Show this message and exit.
    
