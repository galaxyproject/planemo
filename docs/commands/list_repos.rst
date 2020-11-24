
``list_repos`` command
======================================

This section is auto-generated from the help text for the planemo command
``list_repos``. This help message can be generated with ``planemo list_repos
--help``.

**Usage**::

    planemo list_repos [OPTIONS] PROJECT

**Help**

Find all shed repositories in one or more directories and output as yaml.

Currently, a shed repository is considered a directory with a .shed.yml
file.

**Options**::


      --exclude PATH                  Paths to exclude.
      --exclude_from FILE             File of paths to exclude.
      --changed_in_commit_range TEXT  Exclude paths unchanged in git commit range.
      --chunk_count INTEGER           Split output into chunks of this many item and
                                      print --chunk such group.
    
      --chunk INTEGER                 When output is split into --chunk_count
                                      groups, output the group 0-indexedby this
                                      option.
    
      --output TEXT                   File to output to, or - for standard output.
      --help                          Show this message and exit.
    
