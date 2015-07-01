Test Driven Development
=================================

Lets start with a simple wrapper for the BWA_ application (``bwa mem`` in
particular). You can create a new mini-project with a minimal bwa-mem tool
using Planemo's ``project_init`` command.::

    % planemo project_init --template bwa bwa
    % cd bwa

This will create a folder with a ``bwa-mem.xml`` as follows:

.. literalinclude:: writing/bwa-mem_v1.xml
   :language: python
   :emphasize-lines: 8,91-93

.. note:: Highlighted are two relatively recently added enhancements to Galaxy's
    tool XML syntax. The ``check_errors="exit"`` code on the ``command`` block
    will cause Galaxy to use the actual process exit code to determine failure -
    in most cases this is superior to the default Galaxy behavior of checking 
    for the presence of standard error output.

    The ``citations`` block at the bottom will cause Galaxy to generate 
    exportable citations in the tool form and history UIs.

In this form, the tool only accepts a single input. This first thing we will
do to build up this tool is expand that to also allow paired datasets.

Two big ideas behind test-driven development are:

- Write tests first.
- Run the test before you implement the feature. Seeing the initial test failing
  ensures that your feature is actually testing the fix.
