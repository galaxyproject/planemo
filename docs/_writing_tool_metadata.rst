Tool Provided Metadata
======================

Prefer Galaxy datatypes and structured outputs to collect metadata whenever
possible. If an arbitrary number of outputs is needed but no special metadata
must be set, use file-name patterns as described in the
:ref:`dedicated section <multiple-output-files>`.

A metadata manifest is useful when more metadata must be supplied than a file
name pattern can carry, or when you are writing a custom wrapper anyway.

For other cases, a tool can write output metadata to ``galaxy.json`` in the job
working directory. The file name can be changed on the ``outputs`` element:

.. code-block:: xml

    <outputs provided_metadata_file="not_galaxy.json">
        ...
    </outputs>

A tool can generate the file from a config file:

.. code-block:: xml

    <command>cp '$c1' galaxy.json</command>
    <configfiles>
        <configfile name="c1">{"out1": {"name": "Result"}}</configfile>
    </configfiles>

Starting with profile 26.2, Galaxy loads this file only when the tool opts in.
The four opt-ins are an explicit ``provided_metadata_file`` or
``provided_metadata_style`` on ``outputs``, an output with ``format="auto"``,
or a dataset collector with ``discover_via="tool_provided_metadata"`` (including
the equivalent ``from_provided_metadata="true"``). Without one of these, the
file is not loaded and its metadata is silently ignored.

-------------------------------------------
Current Format (Profile 17.09 and Later)
-------------------------------------------

The current format applies to tools with ``profile="17.09"`` or later. It can
also be selected explicitly with ``provided_metadata_style="default"``:

.. code-block:: xml

    <outputs provided_metadata_style="default">
        ...
    </outputs>

The entire metadata file is one JSON object. Each top-level key is the ``name``
of a declared ``data`` or ``collection`` output.

Simple Outputs
--------------

For a simple declared output, the value is its metadata dictionary:

.. code-block:: json

    {
      "out1": {
        "name": "Result",
        "ext": "tabular",
        "info": "Analysis result",
        "dbkey": "hg38",
        "created_from_basename": "input.tsv"
      }
    }

Declared outputs in both the current and legacy formats accept the same fields:

* ``name`` becomes the dataset name.
* ``ext`` sets the dataset format. The output must declare ``format="auto"``
  for this field to take effect.
* ``info`` is the short text shown in the history panel.
* ``dbkey`` sets the genome build.
* ``created_from_basename`` records the basename of the input from which the
  dataset was created.
* ``metadata`` is a dictionary of metadata registered for the datatype.

``auto_format="true"`` sets the format to ``_sniff_`` and does not opt in to
tool-provided metadata. It cannot be combined with ``format``.

See
`tool_provided_metadata_5.xml <https://github.com/galaxyproject/galaxy/blob/c6c3b6df49ff0f17d393fb0c60326d98f3f0251e/test/functional/tools/tool_provided_metadata_5.xml>`__
and
`tool_provided_metadata_9.xml <https://github.com/galaxyproject/galaxy/blob/c6c3b6df49ff0f17d393fb0c60326d98f3f0251e/test/functional/tools/tool_provided_metadata_9.xml>`__.

Discovered Datasets
-------------------

For discovered datasets, the output value contains a ``datasets`` list.
Each entry identifies a discovered file by its basename in ``filename``:

.. code-block:: json

    {
      "report": {
        "datasets": [
          {
            "filename": "sample.report.tsv",
            "name": "Sample report",
            "ext": "tabular",
            "info": "Report for sample",
            "dbkey": "hg38"
          }
        ]
      }
    }

An entry may contain ``filename``, ``designation``, ``name``, ``ext``,
``info``, ``dbkey``, ``metadata``, and ``extra_files``. The ``metadata`` value
is a dictionary of metadata registered for the datatype.

For a discovered dataset, supplying a truthy ``metadata`` dictionary replaces
automatic datatype metadata detection; for a declared output, Galaxy runs
automatic detection and applies the supplied values around it.

Discovery is still declared on the output element:

.. code-block:: xml

    <outputs provided_metadata_style="default">
        <data name="report" format="txt">
            <discover_datasets pattern="(?P&lt;designation&gt;.+)\.report\.tsv"
                               visible="true" />
        </data>
    </outputs>

With pattern discovery, the ``designation`` comes from the regular expression's
named group. To take ``designation`` from each JSON entry instead, use
provided-metadata discovery on a plain ``data`` output:

.. code-block:: xml

    <outputs>
        <data name="report">
            <discover_datasets from_provided_metadata="true" visible="true" />
        </data>
    </outputs>

Examples include
`plain discovery <https://github.com/galaxyproject/galaxy/blob/c6c3b6df49ff0f17d393fb0c60326d98f3f0251e/test/functional/tools/tool_provided_metadata_6.xml>`__,
`provided-metadata discovery <https://github.com/galaxyproject/galaxy/blob/c6c3b6df49ff0f17d393fb0c60326d98f3f0251e/test/functional/tools/tool_provided_metadata_7.xml>`__,
`datatype metadata <https://github.com/galaxyproject/galaxy/blob/c6c3b6df49ff0f17d393fb0c60326d98f3f0251e/test/functional/tools/tool_provided_metadata_10.xml>`__,
and
`extra files <https://github.com/galaxyproject/galaxy/blob/c6c3b6df49ff0f17d393fb0c60326d98f3f0251e/test/functional/tools/tool_provided_metadata_12.xml>`__.

Collection Elements
-------------------

Collection elements can be described with nested ``elements``:

.. code-block:: json

    {
      "list_output": {
        "elements": [
          {
            "name": "oe1",
            "elements": [
              {"name": "ie1", "filename": "oe1_ie1.fq"},
              {"name": "ie2", "filename": "oe1_ie2.fq"}
            ]
          }
        ]
      }
    }

The equivalent flat form uses ``datasets`` and one ``identifier_<level>`` key
for each nesting level:

.. code-block:: json

    {
      "list_output": {
        "datasets": [
          {
            "identifier_0": "oe1",
            "identifier_1": "ie1",
            "filename": "oe1_ie1.fq"
          },
          {
            "identifier_0": "oe1",
            "identifier_1": "ie2",
            "filename": "oe1_ie2.fq"
          }
        ]
      }
    }

Galaxy converts nested ``elements`` to the flat representation internally,
using each element's ``name`` as its level identifier. The collection declares
that its elements come from provided metadata:

.. code-block:: xml

    <outputs>
        <collection name="list_output" type="list:list">
            <discover_datasets from_provided_metadata="true"
                               ext="fastqsanger"
                               visible="true" />
        </collection>
    </outputs>

Do not mix ``datasets`` and ``elements`` in one output entry. Galaxy consults
``elements`` only when ``datasets`` is absent or empty.

See the
`nested elements example <https://github.com/galaxyproject/galaxy/blob/c6c3b6df49ff0f17d393fb0c60326d98f3f0251e/test/functional/tools/collection_creates_dynamic_nested_from_json_elements.xml>`__
and the
`flat datasets example <https://github.com/galaxyproject/galaxy/blob/c6c3b6df49ff0f17d393fb0c60326d98f3f0251e/test/functional/tools/collection_creates_dynamic_nested_from_json.xml>`__.

Output Failures and Reserved Keys
---------------------------------

Setting ``"failed": true`` on any output entry fails the whole job. One entry
carrying the key is enough:

.. code-block:: json

    {
      "out1": {
        "failed": true
      }
    }

``__unnamed_outputs`` is reserved for Galaxy's internal upload handling. Tool
authors should not write it.

-------------------------------------------
Legacy Format (Profiles Before 17.09)
-------------------------------------------

The legacy format applies to profiles older than 17.09. It can be selected
explicitly with ``provided_metadata_style="legacy"``:

.. code-block:: xml

    <outputs provided_metadata_style="legacy">
        ...
    </outputs>

The file uses JSON Lines: one JSON dictionary per line. Every dictionary has a
``type`` key.

Discovered Datasets (Legacy)
----------------------------

A ``new_primary_dataset`` entry describes a discovered dataset. ``filename`` is
the basename of the discovered file:

.. code-block:: json

    {"type": "new_primary_dataset", "filename": "sample.tsv", "name": "Sample", "ext": "tabular", "dbkey": "hg38", "info": "Sample output"}

An entry may also contain ``metadata``, a dictionary of metadata registered for
the datatype, and ``extra_files``.

See
`tool_provided_metadata_2.xml <https://github.com/galaxyproject/galaxy/blob/c6c3b6df49ff0f17d393fb0c60326d98f3f0251e/test/functional/tools/tool_provided_metadata_2.xml>`__,
`tool_provided_metadata_3.xml <https://github.com/galaxyproject/galaxy/blob/c6c3b6df49ff0f17d393fb0c60326d98f3f0251e/test/functional/tools/tool_provided_metadata_3.xml>`__,
and
`tool_provided_metadata_11.xml <https://github.com/galaxyproject/galaxy/blob/c6c3b6df49ff0f17d393fb0c60326d98f3f0251e/test/functional/tools/tool_provided_metadata_11.xml>`__.

Declared Outputs
----------------

A ``dataset`` entry supplies metadata for an explicitly declared output. There
are two ways to identify that output.

``dataset_id`` is the database ID. In practice, tools write
``$out1.dataset.dataset.id`` through a config file:

.. code-block:: json

    {"type": "dataset", "dataset_id": 123, "name": "Result", "ext": "tabular"}

``dataset`` is the output file's basename. Galaxy resolves it by name to that
output's dataset ID:

.. code-block:: json

    {"type": "dataset", "dataset": "${os.path.basename(str($out1))}", "name": "Result", "ext": "tabular"}

These entries accept the declared-output fields described under `Simple Outputs`_.
The ``"failed": true`` behavior described above applies to both
``dataset`` and ``new_primary_dataset`` entries and fails the whole job.

See the ``dataset_id`` examples
`tool_provided_metadata_1.xml <https://github.com/galaxyproject/galaxy/blob/c6c3b6df49ff0f17d393fb0c60326d98f3f0251e/test/functional/tools/tool_provided_metadata_1.xml>`__
and
`tool_provided_metadata_8.xml <https://github.com/galaxyproject/galaxy/blob/c6c3b6df49ff0f17d393fb0c60326d98f3f0251e/test/functional/tools/tool_provided_metadata_8.xml>`__,
and the basename example
`tool_provided_metadata_4.xml <https://github.com/galaxyproject/galaxy/blob/c6c3b6df49ff0f17d393fb0c60326d98f3f0251e/test/functional/tools/tool_provided_metadata_4.xml>`__.