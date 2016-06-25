As shown at the beginning of this section, the command ``seqtk seq`` generates
a help message for the ``seq`` command. ``tool_init`` can take that help message and
stick it right in the generated tool file using the ``help_from_command`` option.

Generally command help messages aren't exactly appropriate for tools
since they mention argument names and simillar details that are abstracted
away by the tool - but they can be an excellent place to start.

The following Planemo's ``tool_init`` call has been enhanced to use ``--help_from_command``.

