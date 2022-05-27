The Basics
====================================================

.. include:: _writing_using_seqtk.rst

For fully featured Seqtk wrappers check out Helena Rasche's
`wrappers <https://github.com/galaxyproject/tools-iuc/tree/master/tools/seqtk>`__
on GitHub.

Galaxy tool files are just XML files, so at this point one could
open a text editor and start writing the tool. Planemo has a command
``tool_init`` to quickly generate some of the boilerplate XML, so let's
start by doing that.

::

    $ planemo tool_init --id 'seqtk_seq' --name 'Convert to FASTA (seqtk)'

The ``tool_init`` command can take various complex arguments - but the two
most basic ones are shown above ``--id`` and ``--name``. Every Galaxy tool
needs an ``id`` (this is a short identifier used by Galaxy itself to identify the
tool) and a ``name`` (this is displayed to the Galaxy user and should be a short
description of the tool). A tool's ``name`` can have whitespace but its ``id``
must not.

The above command will generate the file ``seqtk_seq.xml`` - which looks
like this.

.. literalinclude:: writing/seqtk_seq_v1.xml
   :language: xml

This tool file has the common sections required for a Galaxy tool but you will
still need to open up the editor and fill out the command template, describe
input parameters, tool outputs, write a help section, etc.

The ``tool_init`` command can do a little bit better than this as well. We can
use the test command we tried above ``seqtk seq -A 2.fastq > 2.fasta`` as
an example to generate a command block by specifing the inputs and the outputs
as follows.

::

    $ planemo tool_init --force \
                        --id 'seqtk_seq' \
                        --name 'Convert to FASTA (seqtk)' \
                        --requirement seqtk@1.2 \
                        --example_command 'seqtk seq -A 2.fastq > 2.fasta' \
                        --example_input 2.fastq \
                        --example_output 2.fasta

This will generate the following XML file - which now has correct
definitions for the input and output as well as an actual command template.

.. literalinclude:: writing/seqtk_seq_v2.xml
   :language: xml
   :emphasize-lines: 5-13

.. include:: _writing_from_help_command.rst

::

    $ planemo tool_init --force \
                        --id 'seqtk_seq' \
                        --name 'Convert to FASTA (seqtk)' \
                        --requirement seqtk@1.2 \
                        --example_command 'seqtk seq -A 2.fastq > 2.fasta' \
                        --example_input 2.fastq \
                        --example_output 2.fasta \
                        --test_case \
                        --cite_url 'https://github.com/lh3/seqtk' \
                        --help_from_command 'seqtk seq'

In addition to demonstrating ``--help_from_command``, this demonstrates generating
a test case from our example with ``--test_case`` and adding a citation for the
underlying tool. The resulting tool XML file is:

.. literalinclude:: writing/seqtk_seq_v3.xml
   :language: xml
   :emphasize-lines: 14-55

At this point we have a fairly a functional Galaxy tool with test and help. This was
a pretty simple example - usually you will need to put more work into the tool
to get to this point - ``tool_init`` is really just designed to get you
started.

Now lets lint and test the tool we have developed. The Planemo's ``lint`` (or
just ``l``) command will review tool for XML validity, obvious mistakes, and
compliance with IUC_ best practices.

::

    $ planemo l
    Linting tool /opt/galaxy/tools/seqtk_seq.xml
    Applying linter tests... CHECK
    .. CHECK: 1 test(s) found.
    Applying linter output... CHECK
    .. INFO: 1 outputs found.
    Applying linter inputs... CHECK
    .. INFO: Found 1 input parameters.
    Applying linter help... CHECK
    .. CHECK: Tool contains help section.
    .. CHECK: Help contains valid reStructuredText.
    Applying linter general... CHECK
    .. CHECK: Tool defines a version [0.1.0].
    .. CHECK: Tool defines a name [Convert to FASTA (seqtk)].
    .. CHECK: Tool defines an id [seqtk_seq].
    Applying linter command... CHECK
    .. INFO: Tool contains a command.
    Applying linter citations... CHECK
    .. CHECK: Found 1 likely valid citations.


By default ``lint`` will find all the tools in your current working directory,
but we could have specified a particular tool with ``planemo lint
seqtk_seq.xml``.

Next we can run our tool's functional test with the ``test`` (or just ``t``)
command. This will print a lot of output (as it starts a Galaxy instance) but should
ultimately reveal our one test passed.

.. _DOI: http://www.doi.org/
.. _IUC: https://wiki.galaxyproject.org/IUC
