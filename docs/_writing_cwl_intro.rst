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
                        --example_output 2.fasta
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
    Applying linter docker_image... WARNING
    .. WARNING: Tool does not specify a DockerPull source.
    Applying linter new_draft... CHECK
    .. INFO: Modern CWL version [cwl:draft-3]
    Failed linting

Here the linting failed because we have not yet defined a Docker image for the
the tool. A later revision of this document will cover specifying a Docker image
for this tool with the ``--container`` argument and discuss defining more
parameters for this tool.

For more information on the Common Workflow Language check out the Draft 3
`User Guide`_ and Specification_.

.. _YAML: http://yaml.org/
.. _User Guide: http://www.commonwl.org/draft-3/UserGuide.html
.. _Specification: http://www.commonwl.org/draft-3/CommandLineTool.html

