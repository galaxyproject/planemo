Test-Driven Development
=================================

----------------------------------------------------------------
An Example Tool - BWA
----------------------------------------------------------------

.. note:: This section assumes you have BWA_ installed. If you are using
    a planemo appliance or have Homebrew/linuxbrew installed a different way
    The following command will install BWA_.

    ::

        $ brew install homebrew/science/bwa

Lets start with a simple wrapper for the BWA_ application (``bwa mem`` in
particular). You can create a new mini-project with a minimal bwa-mem tool
using Planemo's ``project_init`` command.::

    $ planemo project_init --template bwa bwa
    $ cd bwa

This will create a folder with a ``bwa-mem.xml`` as follows:

.. literalinclude:: writing/bwa-mem_v1.xml
   :language: xml
   :emphasize-lines: 8,91-93

.. note:: Highlighted are two relatively recently added enhancements to Galaxy's
    tool XML syntax. The ``check_errors="exit"`` code on the ``command`` block
    will cause Galaxy to use the actual process exit code to determine failure -
    in most cases this is superior to the default Galaxy behavior of checking 
    for the presence of standard error output.

    The ``citations`` block at the bottom will cause Galaxy to generate 
    exportable citations in the tool form and history UIs.

----------------------------------------------------------------
Improved Input Handling via Test-Driven Development
----------------------------------------------------------------

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

    $ cd test-data
    $ bwa index -a is bwa-mem-mt-genome.fa
    $ bwa mem bwa-mem-mt-genome.fa bwa-mem-fastq1.fq bwa-mem-fastq2.fq | \
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
        <param name="input_type" value="paired_collection" />
        <output name="bam_output" file="bwa-aln-test2.bam" ftype="bam" lines_diff="2" />
    </test>

We know that we will want to specify the input datasets a paired collection
(see the collections documentation in this document for more information) and
that we will need to have a way to allow the user to specify they are
submitting a paired collection instead of a single input. This is where the
``fastq_input`` and ``input_type`` variables above cam from.

Next run ``planemo lint bwa-mem.xml`` to verify the tool doesn't have any
obvious defects. Once the XML is valid - use planemo ``test`` or just ``t``
for short to verify this new test is failing.

::

    $ planemo t bwa-mem.xml
    ... < bunch of output >
    bwa_mem_test[0]: passed
    bwa_mem_test[1]: failed
    $

Here you can see this second new test is failing - that is good!

.. info: You can open the file ``tool_test_output.html`` in Firefox or your
    favorite web browser to see full details for all tests executed.

        $ firefox tool_test_output.html

The fix is to create a conditional allowing the user to specify an input type.
When modifying the tool and retesting - try passing the ``--failed`` flag to
``planemo test`` - it will speed things up by only rerunning tests have
already failed.

::

    $ planemo t --failed bwa-mem.xml

If you are comfortable with Galaxy tool development - try modifying the tool.

.. info: *Hints*
    * You will need to use a relatively new ``data_collection`` ``param`` type.
      It accepts many of the same attributes as ``data`` parameters (e.g. see
      ``input_fastq1``) but you will need to specify a ``collection_type`` of 
      ``paired``.
    * To access the ``data_collection`` parameter parts in the ``command`` 
      block - use ``$collection_param.forward`` and
      ``$collection_param.reverse``.

Once you get the new test case passing with the ``--failed`` parameter - try
running all the tests again to ensure you didn't break the original test.

::

    $ planemo t bwa-mem.xml
    ... < bunch of output >
    bwa_mem_test[0]: passed
    bwa_mem_test[1]: passed
    $


One possible implementation for tests is as follows (changes highlighted)

.. literalinclude:: writing/bwa-mem_v3.xml
   :language: xml
   :emphasize-lines: 20-25,32-43,58-68

.. note:: *Exercise:* The devteam mappers allow users to specify both a paired
   collection or individual datasets (i.e. using two ``data`` parameters). Extend the above ``conditional`` to allow that. Remember to write your test
   case first and make sure it fails.

   *Hint:* You should not require additional inputs or outputs to do this.

----------------------------------------------------------------
Adding More Parameters
----------------------------------------------------------------

Next up, lets add some of BWA_'s optional parameters to our tool - these
parameters are outlined in the example tool's ``help`` section. To speed this
up and demonstrate another feature of Galaxy - the next test will test the
command-line generated by Galaxy instead of the exact outputs. Not requiring a
complete set of outputs for each test case is convient because it can speed
development and allow testing more parameter combinations but there are
certain tools and certain parameters where exact outputs are impossible to
pre-determine.

Lets start with ``algorithm`` parameter
``-k INT        minimum seed length [19]``. Again, lets do this test first!

::

    <test>
        <param name="fastq_input1" value="bwa-mem-fastq1.fq" />
        <param name="ref_file" value="bwa-mem-mt-genome.fa" />
        <param name="set_algorithm_params" value="true" />
        <param name="k" value="20" />
        <assert_command>
            <has_text text="-k 20" />
        </assert_command>
    </test>

Continuing our pattern - lets ensure this new test fails before implementing
the ``k`` parameter.

::

    $ planemo t bwa-mem.xml
    ... < bunch of output >
    bwa_mem_test[0]: passed
    bwa_mem_test[1]: passed
    bwa_mem_test[2]: failed
    $

Reviewing the output - indeed this new test failed as expected (``did not
contain expected text '-k 20'``). Now let's implement the ``k`` parameter and
use ``planemo t --failed bwa-mem.xml`` to ensure our implementation is
correct.

An example tool with this test and passing.

.. literalinclude:: writing/bwa-mem_v5.xml
   :language: xml
   :emphasize-lines: 18-22,48-56,82-90

The tool also demonstrates the new ``argument`` option on ``param`` tag. These
work a lot like specifying a parameter ``name`` argument - but Galaxy will
describe the underlying application argument in the GUI and API - which
may be helpful for power users and external applications. 

.. note:: *Exercise:* Implement a few more algorithm parameters and start 
    another ``Scoring`` section. Extend the above test case as you go.

.. note:: *Exercise:* Extend the first test case to verify by default none
    of these parameters are present in the command. Use the ``not_has_text``
    tag to do this (e.g. ``<not_has_text text="-k 20">``).

.. note:: A full list of the current assertion elements like these that are 
    allowed can be found `on the tool syntax page
    <https://wiki.galaxyproject.org/Admin/Tools/ToolConfigSyntax#A.3Cassert_contents.3E_tag_set_.28functional_tests.29>`__.

    In additon to the assertion-based testing of the command, the jobs standard
    output and standard error can be checked using ``assert_stdout`` and 
    ``assert_stderr`` respectively - paralleling the ``assert_command`` tag.

    See the example tool `job_properties.xml
    <https://github.com/galaxyproject/galaxy/blob/dev/test/functional/tools/job_properties.xml>`__
    for an example of this.

----------------------------------------------------------------
Publishing to Tool Shed
----------------------------------------------------------------

Publishing to the Tool Shed
====================================================

Now that the tool is working and useful - it is time to publish it to the tool
shed. The `Galaxy Tool Shed`_ (referred to colloquially in Planemo as the
"shed") can store Galaxy tools, dependency definitions, and workflows among
other Galaxy artifacts.

-------------------------------------------------
Configuring a Shed Account
-------------------------------------------------

The `planemo <http://planemo.readthedocs.org/en/latest/appliance.html>`__
appliance comes pre-configured with a local Tool Shed and planemo is
configured to talk to it via ``~/.planemo.yml``. Check out the `publishing docs
<http://planemo.readthedocs.org/en/latest/publishing.html>`__ for information
on setting up this ``~/.planemo.yml`` file on your development environment.
 
-------------------------------------------------
Creating a Repository
-------------------------------------------------

Planemo can be used to used to publish "repositories" to the Tool Shed. A
single GitHub repository or locally managed directory of tools may correspond
to any number of Tool Shed repositories. Planemo maps files to Tool Shed
repositories using a special file called ``.shed.yml``.

From a directory containing tools the ``shed_init``
`command <http://planemo.readthedocs.org/en/latest/commands.html#shed-init-command>`__
can be used to bootstrap a new ``.shed.yml`` file.

::

    planemo shed_init --name=bwa-bwa
                      --owner=planemo
                      --description=bwa-mem
                      --long_description="BWA MEM: Long and medium read mapper"
                      --category="Next Gen Mappers"

There is not a lot of magic happening here, this file could easily be created
directly with a text editor - but the command has a ``--help`` to assist you
and does some very basic validation.

More information on ``.shed.yml`` can be found as part of the IUC `best
practice documentation
<http://galaxy-iuc-standards.readthedocs.org/en/latest/best_practices/shed_yml.html>`__.

After reviewing ``.shed.yml``, this configuration file and relevant shed
artifacts can be quickly linted using the following command.

::

    planemo shed_lint --tools

Once the details the ``.shed.yml`` are set and it is time to create the remote
repository and upload artifacts to it - the following two commands can be used
- the first only needs to be run once and creates the repository based the
metadata in ``.shed.yml`` and the second uploads your actual artifacts to it.

::

    planemo shed_create --shed_target local

Now navigate to the local tool shed (likely at `http://localhost:9009/
<http://localhost:9009/>`__). You can login with login ``planemo@test.com``
and password ``planemo``.

-------------------------------------------------
Updating a Repository
-------------------------------------------------

::

    planemo shed_update --shed_target local

Once tools and reqiured dependency files have been published to the tool shed,
the actual shed dependencies can be automatically and installed and tool 
tests ran using the command::

    planemo shed_serve --shed_target local

Once your artifacts are ready for publication to the main Tool Shed, the
following commands to create a repository there and populate it with your
repository contents.

::

    planemo shed_create

The planemo machine isn't preconfigured to allow publishing to the main tool
shed so this command will not work. See the more complete `publishing docs
<http://planemo.readthedocs.org/en/latest/publishing.html>`__ for full details
about how to setup Planemo to publish to the main and test tool shed - the
process is very similar.
