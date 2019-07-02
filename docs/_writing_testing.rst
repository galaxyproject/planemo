Test-Driven Development
=================================

----------------------------------------------------------------
An Example Tool - BWA
----------------------------------------------------------------

To get started let's install BWA. Here we are going to use ``conda`` to
install BWA - but however you obtain it should be fine.

::

    $ conda install --force -c conda-forge -c bioconda bwa
        ... bwa installation ...
    $ bwa
    Program: bwa (alignment via Burrows-Wheeler transformation)
    Version: 0.7.13-r1126
    Contact: Heng Li <lh3@sanger.ac.uk>

    Usage:   bwa <command> [options]

    Command: index         index sequences in the FASTA format
             mem           BWA-MEM algorithm
             fastmap       identify super-maximal exact matches
             pemerge       merge overlapping paired ends (EXPERIMENTAL)
             aln           gapped/ungapped alignment
             samse         generate alignment (single ended)
             sampe         generate alignment (paired ended)
             bwasw         BWA-SW for long queries

             shm           manage indices in shared memory
             fa2pac        convert FASTA to PAC format
             pac2bwt       generate BWT from PAC
             pac2bwtgen    alternative algorithm for generating BWT
             bwtupdate     update .bwt to the new format
             bwt2sa        generate SA from BWT and Occ

    Note: To use BWA, you need to first index the genome with `bwa index'.
          There are three alignment algorithms in BWA: `mem', `bwasw', and
          `aln/samse/sampe'. If you are not sure which to use, try `bwa mem'
          first. Please `man ./bwa.1' for the manual.


Alternatively you can use Homebrew/linuxbrew to install it:

    ::

        $ brew install homebrew/science/bwa

Lets start with a simple wrapper for the BWA application (``bwa mem`` in
particular). You can create a new mini-project with a minimal bwa-mem tool
using Planemo's ``project_init`` command.

::

    $ planemo project_init --template bwa bwa
    $ cd bwa

This will create a folder with a ``bwa-mem.xml`` as follows:

.. literalinclude:: writing/bwa-mem_v1.xml
   :language: xml
   :emphasize-lines: 8,91-93

Highlighted are two features of Galaxy's tool XML syntax.
The ``detect_errors="exit_code"`` on the ``command`` block
will cause Galaxy to use the actual process exit code to determine failure -
in most cases this is superior to the default Galaxy behavior of checking
for the presence of standard error output.

The ``<citations>`` block at the bottom will cause Galaxy to generate
exportable citations in the tool form and history UIs.

----------------------------------------------------------------
Improved Input Handling via Test-Driven Development
----------------------------------------------------------------

In this form, the tool only accepts a single input. The first thing we will
do is to expand the tool to also allow paired datasets.

.. note:: Two big ideas behind test-driven development are:

  - Write a failing test first.
  - Run the test before you implement the feature. Seeing the initial test failing
    ensures that your feature is actually being tested.

So let's start by generating a test output for the two input files (the bootstrapped
example includes two fastq input files to work with ``bwa-mem-fastq1.fq`` and
``bwa-mem-fastq2.fq``). The following commands will create a bwa index on the
fly, map two input files against it, and build and sort a bam output from the
result - all following the pattern from the command block in the tool.

::

    $ cd test-data
    $ bwa index -a is bwa-mem-mt-genome.fa
    $ bwa mem bwa-mem-mt-genome.fa bwa-mem-fastq1.fq bwa-mem-fastq2.fq | \
      samtools view -Sb - > temporary_bam_file.bam && \
      (samtools sort -f temporary_bam_file.bam bwa-aln-test2.bam || samtools sort -o bwa-aln-test2.bam temporary_bam_file.bam)

.. warning:: In many ways this magic is the hardest part of wrapping Galaxy
   tools and is something this tutorial cannot really teach. The command line
   magic required for each tool is going to be different. Developing a Galaxy
   wrapper requires a lot of knowledge of the underlying applications.

.. note:: Sort appears twice in this odd command because two different
   versions of samtools with conflicting syntaxes may happen
   to be on your machine when running this command. Galaxy manages
   versioned dependencies and so the tool itself does not reflect this
   complexity.

The primary result of this is the file ``test-data/bwa-aln-test2.bam``. We
will  now copy and paste the existing test case to add a new test case that
specifies both fastq inputs as a collection and expects this new output.

.. code-block:: xml

    <test>
        <param name="fastq_input">
            <collection type="paired">
                <element name="forward" value="bwa-mem-fastq1.fq" />
                <element name="reverse" value="bwa-mem-fastq2.fq" />
            </collection>
        </param>
        <param name="ref_file" value="bwa-mem-mt-genome.fa" />
        <param name="input_type" value="paired_collection" />
        <output name="bam_output" file="bwa-aln-test2.bam" ftype="bam" lines_diff="2" />
    </test>

We want to specify the input datasets as a paired collection
(see the collections documentation in this document for more information) and
we need to have a way to allow the user to specify they are
submitting a paired collection instead of a single input. This is where the
``fastq_input`` and ``input_type`` variables above came from.

Next run ``planemo l`` to verify the tool doesn't have any
obvious defects. Once the XML is valid - use ``planemo t``
to verify the new test is failing.

::

    $ planemo t
    ... bunch of output ...
    bwa_mem_test[0]: passed
    bwa_mem_test[1]: failed

.. note:: You can run ``$ firefox tool_test_output.html`` to see
  full details of all executed tests.

Here you can see this second new test is failing - that is good!
The fix is to create a conditional allowing the user to specify an input type.
When modifying the tool and retesting - try passing the ``--failed`` flag to
``planemo t`` - it will speed things up by only rerunning tests that have
already failed.

::

    $ planemo t --failed

If you are comfortable with Galaxy tool development - try modifying the tool
to make the failing test pass.

*Hint:*
  * You will need to use the ``data_collection`` ``param`` type.
    It accepts many of the same attributes as ``data`` parameters (e.g. see
    ``input_fastq1``) but you will need to specify a ``collection_type`` of
    ``paired``.
  * To access the ``data_collection`` parameter parts in the ``command``
    block - use ``$collection_param.forward`` and
    ``$collection_param.reverse``.

Once you get the new test case passing with the ``--failed`` parameter - try
running all the tests again to ensure you didn't break the original test.

::

    $ planemo t
    ... bunch of output ...
    bwa_mem_test[0]: passed
    bwa_mem_test[1]: passed

One possible implementation for tests is as follows (sections with changes highlighted).

.. literalinclude:: writing/bwa-mem_v3.xml
   :language: xml
   :emphasize-lines: 19-26,31-43,58-68

.. note:: **Exercise:** The devteam mappers allow users to specify both a paired
   collection or individual datasets (i.e. using two ``data`` parameters).
   Extend the above ``conditional`` to allow that. Remember to write your test
   case first and make sure it fails.

   *Hint:* You should not require additional inputs or outputs to do this.

----------------------------------------------------------------
Adding More Parameters
----------------------------------------------------------------

Next up, let's add some of BWA's optional parameters to our tool - these
parameters are outlined in the example tool's ``help`` section. To speed this
up and demonstrate another feature of Galaxy - the next test will test the
command-line generated by Galaxy instead of the exact outputs. Not requiring a
complete set of outputs for each test case is convenient because it can speed
development and allows testing more parameter combinations. There are
certain tools and certain parameters where exact outputs are impossible to
pre-determine though.

Lets start with ``algorithm`` parameter``-k INT minimum seed length [19]``.
Again, lets do a test first!

.. code-block:: xml

    <test>
        <param name="fastq_input1" value="bwa-mem-fastq1.fq" />
        <param name="ref_file" value="bwa-mem-mt-genome.fa" />
        <param name="set_algorithm_params" value="true" />
        <param name="k" value="20" />
        <assert_command>
            <has_text text="-k 20" />
        </assert_command>
    </test>

Continuing our pattern - let's ensure this new test fails before implementing
the ``k`` parameter.

::

    $ planemo t
    ... bunch of output ...
    bwa_mem_test[0]: passed
    bwa_mem_test[1]: passed
    bwa_mem_test[2]: failed

Reviewing the output - indeed this new test failed as expected (*did not
contain expected text '-k 20'*). Now let's implement the ``k`` parameter and
use ``planemo t --failed`` to ensure our implementation is correct.

An example tool with this test and passing.

.. literalinclude:: writing/bwa-mem_v5.xml
   :language: xml
   :emphasize-lines: 18-22,48-56,82-90

The tool also demonstrates the new ``argument`` option on ``param`` tag. These
work a lot like specifying a parameter ``name`` argument - but Galaxy will
describe the underlying application argument in the GUI and API - which
may be helpful for power users and external applications.

**Exercise 1:** Implement a few more algorithm parameters and start
another ``Scoring`` section. Extend the above test case as you go.

**Exercise 2:** Extend the first test case to verify by default none
of these parameters are present in the command. Use the ``not_has_text``
tag to do this (e.g. ``<not_has_text text="-k 20">``).

**Exercise 3:** Publish the bwa-mem to the local Tool Shed following
the procedure described in the `tutorial <http://planemo.readthedocs.io/en/latest/writing_appliance.html#publishing-to-the-tool-shed>`__.
(Don't forget to alter the commands from the used ``seqtk`` example to ``bwa-mem``.)

*Hint:*

::

    $ planemo shed_init --name=bwa-bwa \
                        --owner=planemo \
                        --description=bwa-mem \
                        --long_description="BWA MEM: Long and medium read mapper" \
                        --category="Next Gen Mappers"


.. note:: A full list of the current assertion elements like these that are
    allowed can be found `on the tool syntax page
    <https://docs.galaxyproject.org/en/latest/dev/schema.html#tool-tests-test-output-assert-contents>`__.

    In additon to the assertion-based testing of the command, the jobs standard
    output and standard error can be checked using ``assert_stdout`` and
    ``assert_stderr`` respectively - paralleling the ``assert_command`` tag.

    See the sample tool `job_properties.xml
    <https://github.com/galaxyproject/galaxy/blob/dev/test/functional/tools/job_properties.xml>`__
    for an example of this.
