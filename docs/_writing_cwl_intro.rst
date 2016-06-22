The Basics
====================================================

.. include:: _writing_using_seqtk.rst

Common Workflow Language tool files are just simple YAML_ files, so at this point
one could just open a text editor and start implementing the tool. Planemo has a
command ``tool_init`` to quickly generate a skeleton to work from, so let's
start by doing that.

::

    $ planemo tool_init --cwl --id 'seqtk_seq' --name 'Convert to FASTA (seqtk)'

The ``tool_init`` command can take various complex arguments - but three two
most basic ones are shown above ``--cwl``, ``--id`` and ``--name``. The ``--cwl``
flag simply tells Planemo to generate a Common Workflow Language tool. ``--id`` is
a short identifier for this tool and it should be unique across all tools.
``--name`` is a short, human-readable name for the the tool - it corresponds
to the ``label`` attribute in the CWL tool document.

The above command will generate the file ``seqtk_seq.cwl`` - which should look
like this.

.. literalinclude:: writing/seqtk_seq_v1.cwl
   :language: yaml

This tool file has the common fields required for a CWL tool with TODO notes,
but you will still need to open up the editor and fill out the command, describe
input parameters, tool outputs, writeup a help section, etc....

The ``tool_init`` command can do a little bit better than this as well. We can
use the test command we tried above ``seqtk seq -a 2.fastq > 2.fasta`` as
an example to generate a command block by specifing the inputs and the outputs
as follows.

::

    $ planemo tool_init --force \
                        --cwl \
                        --id 'seqtk_seq' \
                        --name 'Convert to FASTA (seqtk)' \
                        --example_command 'seqtk seq -a 2.fastq > 2.fasta' \
                        --example_input 2.fastq \
                        --example_output 2.fasta

This will generate the following CWL tool definition - which now has correct
definitions for the input, output, and command specified. These represent a best
guess by planemo, and in most cases will need to be tweaked manually after the
tool is generated.

.. literalinclude:: writing/seqtk_seq_v2.cwl
   :language: yaml

.. include:: _writing_from_help_command.rst

::

    $ planemo tool_init --force \
                        --cwl \
                        --id 'seqtk_seq' \
                        --name 'Convert to FASTA (seqtk)' \
                        --example_command 'seqtk seq -a 2.fastq > 2.fasta' \
                        --example_input 2.fastq \
                        --example_output 2.fasta \
                        --container 'dukegcb/seqtk' \
                        --test_case \
                        --help_from_command 'seqtk seq'

This command generates the following CWL YAML file.

.. literalinclude:: writing/seqtk_seq_v3.cwl
   :language: yaml

.. include:: _writing_lint_intro.rst

::

    $ planemo l seqtk_seq.cwl
    Linting tool /home/john/workspace/planemo/docs/writing/seqtk_seq.cwl
    Applying linter general... CHECK
    .. CHECK: Tool defines a version [0.0.1].
    .. CHECK: Tool defines a name [Convert to FASTA (seqtk)].
    .. CHECK: Tool defines an id [seqtk_seq_v3].
    Applying linter cwl_validation... CHECK
    .. INFO: CWL appears to be valid.
    Applying linter docker_image... CHECK
    .. INFO: Tool will run in Docker image [dukegcb/seqtk].
    Applying linter new_draft... CHECK
    .. INFO: Modern CWL version [cwl:draft-3]

In addition to the actual tool file, a test file will be generated
using the example command and provided test data. The file contents are as
follows:

.. literalinclude:: writing/seqtk_seq_v3_tests.yml
   :language: yaml

This file is a planemo-specific artifact. This file may contain 1 or more
tests - each test is an element of the top-level list. ``tool_init`` will use
the example command to build just one test.

Each test consists of a few parts:

- ``doc`` - this attribute provides a short description for the test.
- ``job`` - this can be the path to a CWL job description or a job
  description embedded right in the test (``tool_init`` builds the latter).
- ``outputs`` - this section describes the expected output for a test. Each
  output ID of the tool or workflow under test can appear as a key. The
  example above just describes expected specific output file contents exactly
  but many more expectations can be described.

The tests described in this file can be run using the planemo ``test`` (or
simply ``t``) command on the original file. By default, planemo will run tool
tests with Galaxy but we can also specify the use of ``cwltool`` (the
reference implementation of CWL) which will be quicker and more robust until
while Galaxy support for the CWL is still in development.

    $ planemo test --no-container --engine cwltool seqtk_seq.cwl
    Enable beta testing mode to test artifact that isn't a Galaxy tool.
    All 1 test(s) executed passed.
    seqtk_seq_0: passed

We can also open up the Galaxy web inteface with this tool loaded
using the ``serve`` (or just ``s``) command.

::

    $ planemo s --cwl seqtk_seq.cwl
    ...
    serving on http://127.0.0.1:9090

Open up http://127.0.0.1:9090 in a web browser to view your new
tool.

For more information on the Common Workflow Language check out the Draft 3
`User Guide`_ and Specification_.

.. _YAML: http://yaml.org/
.. _User Guide: http://www.commonwl.org/draft-3/UserGuide.html
.. _Specification: http://www.commonwl.org/draft-3/CommandLineTool.html

