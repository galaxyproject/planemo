Test Driven Development
=================================

.. note:: This section assumes you have BWA_ installed. If you are using
    a planemo appliance or have Homebrew/linuxbrew installed a different way
    The following command will install BWA_.

    ::

        % brew install homebrew/science/bwa

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

- Write a failing test first.
- Run the test before you implement the feature. Seeing the initial test failing
  ensures that your feature is actually being tested.

So lets start by generating a test output for two input files (the bootstrapped
example includes two fastq input files to work with ``bwa-mem-fastq1.fq`` and
``bwa-mem-fastq2.fq``). The following command will create a bwa index on the 
fly, map two input files against it, and build and sort a bam output from the 
result - all following the pattern from the command block in the tool.

::

    % cd test-data
    % bwa index -a is bwa-mem-mt-genome.fa
    % bwa mem bwa-mem-mt-genome.fa bwa-mem-fastq1.fq bwa-mem-fastq2.fq | \
      samtools view -Sb - > temporary_bam_file.bam && \
      samtools sort -f temporary_bam_file.bam bwa-aln-test2.bam

.. warning: In many ways this magic is the hardest part of wrapping Galaxy 
   tools and is something this tutorial cannot really teach. The command line 
   magic required for each tool is going to be different. Developing a Galaxy
   wrapper requires a lot of knowledge of the underlying application.s

The primary result of this is the file ``test-data/bwa-aln-test2.bam``. We
will  now copy and paste the existing test case to add a new test case that
specifies both fastq inputs as a collection and expects this new output.

::

    <test>
        <param name="fastq_input">
            <collection type="paired">
                <element name="forward" value="bwa-mem-fastq1.fq" />
                <element name="reverse" value="bwa-mem-fastq2.fq" />
            </collection>
        </param>
        <param name="ref_file" value="bwa-mem-mt-genome.fa" />
        <output name="bam_output" file="bwa-aln-test2.bam" ftype="bam" lines_diff="2" />
    </test>

TODO: SHOW TEST FAILING

TODO: SHOW FIX

TODO: RERUN FAILED TEST


.. note:: *Exercise:* The devteam mappers allow users to specify both a paired
   collection or individual datasets (i.e. using two ``data`` parameters). Extend the above ``conditional`` to allow that. Remember to write your test
   case first and make sure it fails.

   *Hint:* You should not require additional inputs or outputs to do this.

TODO: TRANSITION TO PARAMETERS.

TODO: ADD TEST DEMONSTRATING COMMAND-LINE CHECKING

TODO: EXERCISE - MORE PARAMETERS

TODO: OUTLINE TEST ASSERTIONS.

TODO: DESCRIBE OTHER INEXACT TESTING OPTIONS

TODO: DESCRIBE TESTING METADATA, STANDARD OUTPUT, STANDARD ERROR
