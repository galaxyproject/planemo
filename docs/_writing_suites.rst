Macros
======

If your desire is to write a tool for a single relatively simple application
or script - this section should be skipped. If you hope to maintain a
collection of related tools - experience suggests you will realize there is a
lot of duplicated XML to do this well. Galaxy tool XML macros_ can help reduce
this duplication.

Planemo's ``tool_init`` command can be used to generate a macro file
appropriate for suites of tools by using the ``--macros`` flag. Consider the
following variant of the previous ``tool_init`` command (the only difference
is now we are adding the ``--macros`` flag).

::

    $ planemo tool_init --force \
                        --macros \
                        --id 'seqtk_seq' \
                        --name 'Convert to FASTA (seqtk)' \
                        --requirement seqtk@1.2 \
                        --example_command 'seqtk seq -A 2.fastq > 2.fasta' \
                        --example_input 2.fastq \
                        --example_output 2.fasta \
                        --test_case \
                        --help_from_command 'seqtk seq'

This will produce the two files in your current directory instead of just one
- ``seqtk_seq.xml`` and ``macros.xml``.

.. literalinclude:: writing/seqtk_seq_with_macros.xml
   :language: xml
   :emphasize-lines: 2-6,46

.. literalinclude:: writing/seqtk_macros.xml
   :language: xml

As you can see in the above code macros are reusable chunks of XML that make it easier
to avoid duplication and keep your XML concise.

Further reading:

- `Macros syntax <https://wiki.galaxyproject.org/Admin/Tools/ToolConfigSyntax#Reusing_Repeated_Configuration_Elements>`__ on the Galaxy Wiki.
- `GATK tools <https://github.com/galaxyproject/tools-iuc/tree/master/tools/gatk2>`__ (example tools making extensive use of macros)
- `gemini tools <https://github.com/galaxyproject/tools-iuc/tree/master/tools/gemini>`__ (example tools making extensive use of macros)
- `bedtools tools <https://github.com/galaxyproject/tools-iuc/tree/master/tools/bedtools>`__ (example tools making extensive use of macros)
- Macros implementation details - `Pull Request #129 <https://bitbucket.org/galaxy/galaxy-central/pull-request/129/implement-macro-engine-to-reduce-tool/diff>`__ and `Pull Request #140 <https://bitbucket.org/galaxy/galaxy-central/pull-request/140/improvements-to-tool-xml-macroing-system/diff>`__

.. _macros: https://wiki.galaxyproject.org/Admin/Tools/ToolConfigSyntax#Reusing_Repeated_Configuration_Elements

