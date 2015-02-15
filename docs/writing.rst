======================
Building a Galaxy Tool
======================

This guide is going to demontrate building up Galaxy tools wrappers for
commands from Heng Li's Seqtk_ package - a package for processing sequence
data in FASTA and FASTQ files. For fully worked through Seqtk wrappers -
checkout Eric Rasche's `wrappers <https://github.com/galaxyproject/tools-
iuc/tree/master/tools/seqtk>`_ on Github.

To get started lets install Seqtk, download an example FASTQ file, and test
out the a simple Seqtk command - ``seq`` which converts FASTQ files into
FASTA. Here we are going to use ``brew`` to install seqtk - but however you
obtain it should be fine.

::

    % brew tap homebrew/science
    % brew install seqtk
    ==> Installing seqtk from homebrew/homebrew-science
    ==> Downloading https://github.com/lh3/seqtk/archive/73866e7.tar.gz
    ######################################################################## 100.0%
    ==> make
    /home/john/.linuxbrew/Cellar/seqtk/1.0-r68: 3 files, 208K, built in 2 seconds
    % wget https://raw.githubusercontent.com/galaxyproject/galaxy-test-data/master/2.fastq
    % seqtk seq

        Usage:   seqtk seq [options] <in.fq>|<in.fa>

        Options: -q INT    mask bases with quality lower than INT [0]
                 -X INT    mask bases with quality higher than INT [255]
                 -n CHAR   masked bases converted to CHAR; 0 for lowercase [0]
                 -l INT    number of residues per line; 0 for 2^32-1 [0]
                 -Q INT    quality shift: ASCII-INT gives base quality [33]
                 -s INT    random seed (effective with -f) [11]
                 -f FLOAT  sample FLOAT fraction of sequences [1]
                 -M FILE   mask regions in BED or name list FILE [null]
                 -L INT    drop sequences with length shorter than INT [0]
                 -c        mask complement region (effective with -M)
                 -r        reverse complement
                 -A        force FASTA output (discard quality)
                 -C        drop comments at the header lines
                 -N        drop sequences containing ambiguous bases
                 -1        output the 2n-1 reads only
                 -2        output the 2n reads only
                 -V        shift quality by '(-Q) - 33'
    % seqtk seq -a 2.fastq > 2.fasta
    % cat 2.fasta
    >EAS54_6_R1_2_1_413_324
    CCCTTCTTGTCTTCAGCGTTTCTCC
    >EAS54_6_R1_2_1_540_792
    TTGGCAGGCCAAGGCCGATGGATCA
    >EAS54_6_R1_2_1_443_348
    GTTGCTTCTGGCGTGGGTGGGGGGG

Galaxy tool files are just simple XML files, so at this point one could just
open a text editor and start implementing the tool. Planemo has a command
``tool_init`` to quickly generate some of the boilerplate XML, so lets
start by doing that.

::

    % planemo tool_init --id 'seqtk_seq' --name 'Convert to FASTA (seqtk)'

The ``tool_init`` command can take various complex arguments - but the two
most basic ones are shown above ``--id`` and ``--name``. Every Galaxy tool
needs an ``id` (this a short identifier used by Galaxy itself to identify the
tool) and a ``name`` (this is display to the Galaxy user and should be a short
description of the tool). In general ``name`` can have whitespace and ``id``
should not.

The above command will generate the file ``seqtk_seq.xml`` - which should look
like this.

.. literalinclude:: writing/seqtk_seq_v1.xml
   :language: xml

This tool file has the common sections required for Galaxy tool but you will
still need to open up the editor and fill out the command template, describe
input parameters, tool outputs, writeup a help section, etc....

The ``tool_init`` command can do a little bit better than this as well. We can
use the test command we generated above ``seqtk seq -a 2.fastq > 2.fasta`` as
an example to generate a command block by specifing the inputs and the outputs
as follows.

::

    % planemo tool_init --force \
                        --id 'seqtk_seq' \
                        --name 'Convert to FASTA (seqtk)' \
                        --requirement seqtk@1.0-r68 \
                        --example_command 'seqtk seq -a 2.fastq > 2.fasta' \
                        --example_input 2.fastq \
                        --example_output 2.fasta

This will generate the following tool XML file - which now has correct
definitions for the input and output as well as an actual command template.

.. literalinclude:: writing/seqtk_seq_v2.xml
   :language: xml

As shown above the command ``seqtk seq`` generates a help message for this the
``seq`` command. ``tool_init`` can take that help message and stick it right
in the generated tool file using the ``help_from_command`` option. Generally
command help messages aren't exactly appropriate for Galaxy tool wrappers
since they mention argument names and simillar details that are abstracted
away by the tool - but they can be a good place to start.

::

    % planemo tool_init --force \
                        --id 'seqtk_seq' \
                        --name 'Convert to FASTA (seqtk)' \
                        --requirement seqtk@1.0-r68 \
                        --example_command 'seqtk seq -a 2.fastq > 2.fasta' \
                        --example_input 2.fastq \
                        --example_output 2.fasta \
                        --test_case \
                        --help_from_command 'seqtk seq'

.. literalinclude:: writing/seqtk_seq_v3.xml
   :language: xml

At this point we have a fairly a functional tool with test and help. This was
a pretty simple example - usually you will need to put more work into the tool
XML to get to this point - ``tool_init`` is really just designed to get you
started.

Now lets lint and test the tool we have developed. The planemo ``lint`` (or
just ``l``) command will reviews tools for obvious mistakes and compliance
with best practices.

::

    % planemo l
    Linting tool /home/john/test/seqtk_seq.xml
    Applying linter lint_top_level... CHECK
    .. CHECK: Tool defines a version.
    .. CHECK: Tool defines a name.
    .. CHECK: Tool defines an id name.
    Applying linter lint_tests... CHECK
    .. CHECK: 1 test(s) found.
    Applying linter lint_output... CHECK
    .. INFO: 1 output datasets found.
    Applying linter lint_inputs... CHECK
    .. INFO: Found 1 input parameters.
    Applying linter lint_help... CHECK
    .. CHECK: Tool contains help section.
    .. CHECK: Help contains valid reStructuredText.
    Applying linter lint_command... CHECK
    .. INFO: Tool contains a command.
    Applying linter lint_citations... WARNING
    .. WARNING: No citations found, consider adding citations to your tool.

By default ``lint`` will find all the tools in your current working directory,
but we could have specified a particular tool with ``planemo lint
seqtk_seq.xml``. The only warning we received here is telling us the tool
lacks citations. Seqtk_ is unpublished, but if there were a paper to cite we
could have done so by passing in the DOI_ (e.g. ``--doi '10.1101/010538'``).

Next we can run our tool's functional test with the ``test`` (or just ``t``)
command. This will print a lot of output but should ultimately reveal our one
test passed.

::

    % planemo --galaxy_root=/path/to/galaxy t
    ...
    All 1 test(s) executed passed.
    seqtk_seq[0]: passed

Now we can open Galaxy 

::

    % planemo --galaxy_root=/path/to/galaxy s
    ...
    serving on http://127.0.0.1:9090

Open up http://127.0.0.1:9090 in a web browser to view your new tool.

More information:

 * `Galaxy's Tool XML Syntax <https://wiki.galaxyproject.org/Admin/Tools/ToolConfigSyntax>`_
 * `Big List of Tool Development Resources <https://wiki.galaxyproject.org/Develop/ResourcesTools>`_

.. _DOI: http://www.doi.org/
.. _Seqtk: https://github.com/lh3/seqtk
