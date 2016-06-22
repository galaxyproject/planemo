Wrapping a Script
====================================================

Many common bioinformatics applications are available on the `Tool Shed`_
already and so a common development task is to integrate scripts of
various complexity into Galaxy.

Consider the following small Perl script.

.. literalinclude:: writing/gc_content.pl
   :language: perl

One can build a tool for this script as follows and place the script in
the same directory as the tool XML file itself. The special value
``$__tool_directory__`` here refers to the directory your tool lives in.

.. literalinclude:: writing/gc_content.xml
   :language: xml

.. _Tool Shed: http://toolshed.g2.bx.psu.edu/
