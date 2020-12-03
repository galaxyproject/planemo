
``ci_find_tools`` command
======================================

This section is auto-generated from the help text for the planemo command
``ci_find_tools``. This help message can be generated with ``planemo ci_find_tools
--help``.

**Usage**::

    planemo ci_find_tools [OPTIONS] PROJECT

**Help**

Find all tools in one or more directories.

Tools can be chunked up, filtered, etc... to build lists of tools to perform
operations over for continuous integration operations.

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
      --group_tools                   Group tools of the same repository on a single
                                      line.
    
      --help                          Show this message and exit.
    
