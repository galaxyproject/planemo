Tool Provided Metadata
==============================

This stub of a section provides some initial documentation on tool provided
metadata. Galaxy allows datasets to be discovered after a tool has been
executed and allows tools to specify metadata for these datasets. Whenever
possible, Galaxy's datatypes and more structured outputs should
be used to collect metadata.

If an arbitrary number of outputs is needed but no special metadata must be set,
file name patterns can be used to allow Galaxy to discover these datasets.
More information on this can be found in the :ref:`dedicated section <multiple-output-files>`.

The file name patterns described in the above link are nice because they don't
need special instrumenting in the tool wrapper to adapt to Galaxy in general and
can adapt to many existing application's output. When more metadata must be 
supplied or when implementing a custom tool wrapper anyway - it may be beneficial
to build a manifest file.

A tool may also produce a file called ``galaxy.json`` during execution. If 
upon a job's completion this file is populated, Galaxy will expect to find metadata
about outputs in it.

The format of this file is a bit quirky - each line of this file should be a JSON
dictionary. Each such dictionary should contain a ``type`` attribute - this type
may be either ``new_primary_dataset`` or ``dataset``. 

If the ``type`` is ``new_primary_dataset``, the dictionary should contain a 
``filename`` entry with a path to a "discovered dataset". In this case the 
dictionary may contain any of the following entries ``name``, ``dbkey``, ``info``, ``ext``, ``metadata``.

- ``name`` will be used as the output dataset's name
- ``ext`` allows specification of the format of the output (e.g. ``txt``, ``tabular``, ``fastqsanger``, etc...)
- ``dbkey`` allows specifying a genome build for the discovered dataset
- ``info`` is a short text description for each dataset that appears in the history panel
- ``metadata`` this should be a dictionary of key-value pairs for metadata registered with the datatype for this output

Examples of tools using ``new_primary_dataset`` entries:

- `tool_provided_metadata_2.xml <https://github.com/jmchilton/galaxy/blob/2909e74642180bd818019ebdcb62e62f12e56e69/test/functional/tools/tool_provided_metadata_2.xml>`__ demonstrating using the simpler attributes described here.
- `tool_provided_metadata_3.xml <https://github.com/jmchilton/galaxy/blob/2909e74642180bd818019ebdcb62e62f12e56e69/test/functional/tools/tool_provided_metadata_3.xml>`__ demonstrates overridding datatype specified metadata.

The ``type`` of an entry may also be ``dataset``. In this case the metadata 
descriptions describe an explicit output (one with its own corresponding ``output``
element definition). In this case, an entry called ``dataset`` should appear in
the dictionary (in lieu of ``filename`` above) and should be the database id of the 
output dataset. Such entries may contain all of the other fields described above except
``metadata``.

Example tool using a ``dataset`` entry:

- `tool_provided_metadata_1.xml <https://github.com/jmchilton/galaxy/blob/2909e74642180bd818019ebdcb62e62f12e56e69/test/functional/tools/tool_provided_metadata_1.xml>`__
