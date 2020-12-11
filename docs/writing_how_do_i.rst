====================================================
How do I...
====================================================

This section contains a number of smaller topics with links and examples meant
to provide relatively concrete answers for specific tool development
scenarios.

------------------------------------------
\.\.\. deal with index/reference data?
------------------------------------------

Galaxy's concept of `data tables
<https://wiki.galaxyproject.org/Admin/Tools/Data%20Tables>`__ are meant to
provide tools with access reference datasets or index data not tied to
particular histories or users. A common example would be FASTA files for
various genomes or mapper-specific indices of those files (e.g. a BWA index
for the hg19 genome).

Galaxy `data managers
<https://wiki.galaxyproject.org/Admin/Tools/DataManagers>`__ are specialized
tools designed to populate tool data tables.


------------------------------------------
\.\.\. cite tools without an obvious DOI?
------------------------------------------

In the absence of an obvious DOI_, tools may contain embedded BibTeX_ directly.

Futher reading:

- `bibtex.xml <https://github.com/jmchilton/galaxy/blob/dev/test/functional/tools/bibtex.xml>`__ (test tool with a bunch of random examples)
- `bwa-mem.xml <https://github.com/jmchilton/bwa-mem/commit/0425264039950bfd9ded06997a08cc8b4ee1ad8f>`__ (BWA-MEM tool by Anton Nekrutenko demonstrating citation of an arXiv article)
- `macros.xml <https://github.com/galaxyproject/tools-devteam/blob/master/tool_collections/vcflib/macros.xml#L15>`__ (Macros for vcflib tool demonstrating citing a github repository)

--------------------------------------------------
\.\.\. declare a Docker container for my tool?
--------------------------------------------------

Galaxy tools can be decorated to with ``container`` tags indicated Docker
container ids that the tools can run inside of.

The longer term plan for the Tool Shed ecosystem is to be able to
automatically build Docker containers for tool dependency descriptions and
thereby obtain this Docker functionality for free and in a way that is
completely backward compatible with non-Docker deployments.

Further reading:

- `Complete tutorial <https://github.com/apetkau/galaxy-hackathon-2014>`__
  on Github by Aaron Petkau. Covers installing Docker, building a Dockerfile_, publishing to `Docker Hub`_, annotating tools and configuring Galaxy.
- `Another tutorial <https://www.e-biogenouest.org/groups/guggo>`__
  from the Galaxy User Group Grand Ouest.
- Landing page on the `Galaxy Wiki <https://wiki.galaxyproject.org/Admin/Tools/Docker>`__
- Impementation details on `Pull Request #401 <https://bitbucket.org/galaxy/galaxy-central/pull-request/401/allow-tools-and-deployers-to-specify>`__

--------------------------------------------------
\.\.\. do extra validation of parameters?
--------------------------------------------------

Tool parameters support a ``validator`` element (`syntax
<https://wiki.galaxyproject.org/Admin/Tools/ToolConfigSyntax#A.3Cvalidator.3E_tag_set>`__)
to perform validation of a single parameter. More complex validation across
parameters can be performed using arbitrary Python functions using the
``code`` file syntax but this feature should be used sparingly.

Further reading:

- `validator <https://wiki.galaxyproject.org/Admin/Tools/ToolConfigSyntax#A.3Cvalidator.3E_tag_set>`__
  XML tag syntax on the Galaxy wiki.
- `fastq_filter.xml <https://github.com/galaxyproject/tools-devteam/blob/master/tool_collections/galaxy_sequence_utils/fastq_filter/fastq_filter.xml>`__
  (a FASTQ filtering tool demonstrating validator constructs)
- `gffread.xml <https://github.com/galaxyproject/tools-devteam/blob/master/tool_collections/cufflinks/gffread/gffread.xml>`__
  (a tool by Jim Johnson demonstrating using regular expressions with ``validator`` tags)
- `code_file.xml <https://github.com/galaxyproject/galaxy/blob/dev/test/functional/tools/code_file.xml>`__,
  `code_file.py <https://github.com/galaxyproject/galaxy/blob/dev/test/functional/tools/code_file.py>`__
  (test files demonstrating defining a simple constraint in Python across
  two parameters)
- `deseq2 tool <https://github.com/bgruening/galaxytools/tree/master/tools/deseq2>`__
  by Björn Grüning demonstrating advanced ``code`` file validation.

-------------------------------------------------
\.\.\. check input type in command blocks?
-------------------------------------------------

Input data parameters may specify multiple formats. For example

.. code-block:: xml

    <param name="input" type="data" format="fastq,fasta" label="Input" />

If the command-line under construction doesn't require changes based
on the input type - this may just be referenced as ``$input``. However, if the
command-line under construction uses different argument names depending on
type for instance - it becomes important to dispatch on the underlying type.

In this example ``$input.ext`` - would return the short code for the actual
datatype of the input supplied - for instance the string ``fasta`` or
``fastqsanger`` would be valid responses for inputs to this parameter for the
above definition.

While ``.ext`` may sometimes be useful - there are many cases where it is
inappropriate because of subtypes - checking if ``.ext`` is equal to ``fastq``
in the above example would not catch ``fastqsanger`` inputs for instance. To
check if an input matches a type or any subtype thereof - the ``is_of_type``
method can be used. For instance

::

    $input.is_of_type('fastq')

would check if the input is of type ``fastq`` or any derivative types such as
``fastqsanger``.

- `Pull Request 457 <https://bitbucket.org/galaxy/galaxy-central/pull-request/457/allow-cheetah-tool-templates-to-reason/diff>`__

-------------------------------------------------
\.\.\. handle arbitrary output data formats?
-------------------------------------------------

If the output format of a tool's output cannot be known ahead of time,
Galaxy can be instructed to "sniff" the output and determine the data type
using the same method used for uploads. Adding the ``auto_format="true"``
attribute to a tool's output enables this.

.. code-block:: xml

    <output name="out1" auto_format="true" label="Auto Output" />

- `output_auto_format.xml <https://github.com/galaxyproject/galaxy/blob/dev/test/functional/tools/output_auto_format.xml>`__

-------------------------------------------------
\.\.\. determine the user submitting a job?
-------------------------------------------------

The variable ``$__user_email__`` (as well as ``$__user_name__`` and
``$__user_id__``) is available when building up your command in
the tool's ``<command>`` block. The following tool demonstrates the use of
this and a few other special parameters available to all tools.

- `special_params.xml <https://github.com/galaxyproject/galaxy/blob/dev/test/functional/tools/special_params.xml>`__

------------------------------------------
\.\.\. test with multiple value inputs?
------------------------------------------

To write tests that supply multiple values to a ``multiple="true"`` ``select`` or ``data`` parameter - simply specify the multiple values as a comma seperated list.

Here are examples of each:

- `multi_data_param.xml <https://github.com/galaxyproject/galaxy/blob/dev/test/functional/tools/multi_data_param.xml>`__
- `muti_select.xml <https://github.com/galaxyproject/galaxy/blob/dev/test/functional/tools/multi_select.xml>`__

------------------------------------------
\.\.\. test dataset collections?
------------------------------------------

Here are some examples of testing tools that consume collections with ``type="data_collection"`` parameters.

- `collection_paired_test.xml <https://github.com/galaxyproject/galaxy/blob/dev/test/functional/tools/collection_paired_test.xml>`__
- `collection_mixed_param.xml <https://github.com/galaxyproject/galaxy/blob/dev/test/functional/tools/collection_mixed_param.xml>`__
- `collection_nested_param.xml <https://github.com/galaxyproject/galaxy/blob/dev/test/functional/tools/collection_nested_test.xml>`__

Here are some examples of testing tools that produce collections with ``output_collection`` elements.

- `collection_creates_list.xml <https://github.com/galaxyproject/galaxy/blob/dev/test/functional/tools/collection_creates_list.xml>`__
- `collection_creates_list_2.xml <https://github.com/galaxyproject/galaxy/blob/dev/test/functional/tools/collection_creates_list_2.xml>`__
- `collection_creates_pair.xml <https://github.com/galaxyproject/galaxy/blob/dev/test/functional/tools/collection_creates_pair.xml>`__
- `collection_creates_pair_from_type.xml <https://github.com/galaxyproject/galaxy/blob/dev/test/functional/tools/collection_creates_pair_from_type.xml>`__

------------------------------------------
\.\.\. test discovered datasets?
------------------------------------------

Tools which dynamically `discover datasets
<https://wiki.galaxyproject.org/Admin/Tools/Multiple%20Output%20Files#Number_of_Output_datasets_cannot_be_determined_until_tool_run>`__
after the job is complete, either using the ``<discovered_datasets>`` element,
the older default pattern approach (e.g. finding files with names like
``primary_DATASET_ID_sample1_true_bam_hg18``), or the undocumented
``galaxy.json`` approach can be tested by placing ``discovered_dataset``
elements beneath the corresponding ``output`` element with the ``designation``
corresponding to the file to test.

.. code-block:: xml

    <test>
      <param name="input" value="7" />
      <output name="report" file="example_output.html">
        <discovered_dataset designation="world1" file="world1.txt" />
        <discovered_dataset designation="world2">
          <assert_contents>
            <has_line line="World Contents" />
          </assert_contents>
        </discovered_dataset>
      </output>
    </test>

The test examples distributed with Galaxy demonstrating dynamic discovery and
the testing thereof include:

- `multi_output.xml <https://github.com/galaxyproject/galaxy/blob/dev/test/functional/tools/multi_output.xml>`__
- `multi_output_assign_primary.xml <https://github.com/galaxyproject/galaxy/blob/dev/test/functional/tools/multi_output_assign_primary.xml>`__
- `multi_output_configured.xml <https://github.com/galaxyproject/galaxy/blob/dev/test/functional/tools/multi_output_configured.xml>`__

------------------------------------------
\.\.\. test composite dataset contents?
------------------------------------------

Tools which consume Galaxy `composite datatypes
<https://wiki.galaxyproject.org/Admin/Datatypes/Composite%20Datatypes>`__ can
generate test inputs using the ``composite_data`` element demonstrated by the
following tool.

- `composite.xml <https://github.com/galaxyproject/galaxy/blob/dev/test/functional/tools/composite.xml>`__

Tools which produce Galaxy `composite datatypes
<https://wiki.galaxyproject.org/Admin/Datatypes/Composite%20Datatypes>`__ can
specify tests for the individual output files using the ``extra_files`` element
demonstrated by the following tool.

- `composite_output.xml <https://github.com/galaxyproject/galaxy/blob/dev/test/functional/tools/composite_output.xml>`__
- `macs_wrapper.xml <https://github.com/galaxyproject/tools-devteam/blob/master/tools/macs/macs_wrapper.xml>`__

------------------------------------------
\.\.\. test index (\.loc) data?
------------------------------------------

There is an idiom to supply test data for index during tests using Planemo_.

To create this kind of test, one needs to provide a
``tool_data_table_conf.xml.test`` beside your tool's
``tool_data_table_conf.xml.sample`` file that specifies paths to test ``.loc``
files which in turn define paths to the test index data. Both the ``.loc``
files and the ``tool_data_table_conf.xml.test`` can use the value
``${__HERE__}`` which will be replaced with the path to the directory the file
lives in. This allows using relative-like paths in these files which is needed
for portable tests.

An example commit demonstrating the application of this approach to a Picard_
tool can be found `here <https://github.com/jmchilton/picard/commit/4df8974384081ee1bb0f97e1bb8d7f935ba09d73>`__.

These tests can then be run with the Planemo `test command
<http://planemo.readthedocs.org/en/latest/commands.html#test-command>`__.


------------------------------------------
\.\.\. test exit codes?
------------------------------------------

A ``test`` element can check the exit code of the underlying job using the
``check_exit_code="n"`` attribute.

- `job_properties.xml <https://github.com/galaxyproject/galaxy/blob/dev/test/functional/tools/job_properties.xml>`__

------------------------------------------
\.\.\. test failure states?
------------------------------------------

Normally, all tool test cases described by a ``test`` element are expected to
pass - but on can assert a job should fail by adding ``expect_failure="true"``
to the ``test`` element.

- `job_properties.xml <https://github.com/galaxyproject/galaxy/blob/dev/test/functional/tools/job_properties.xml>`__

------------------------------------------
\.\.\. test output filters work?
------------------------------------------

If your tool contains ``filter`` elements, you can't verify properties of outputs
that are filtered out and do not exist. The ``test`` element may contain an
``expect_num_outputs`` attribute to specify the expected number of outputs, this
can be used to verify that outputs not listed are expected to be filtered out during
tool execution.

- `output_filter.xml <https://github.com/galaxyproject/galaxy/blob/dev/test/functional/tools/output_filter.xml>`__

------------------------------------------
\.\.\. test metadata?
------------------------------------------

Output metadata can be checked using ``metadata`` elements in the XML
description of the ``output``.

- `metadata.xml <https://github.com/galaxyproject/galaxy/blob/dev/test/functional/tools/metadata.xml>`__

--------------------------------------------------------------------
\.\.\. test tools installed in an existing Galaxy instance?
--------------------------------------------------------------------

Do not use planemo, Galaxy should be used to test its tools directly.
The following two commands can be used to test Galaxy tools in an existing
instance.

::

    $ sh run_tests.sh --report_file tool_tests_shed.html --installed


This above command specifies the ``--installed`` flag when calling
``run_tests.sh``, this tells the test framework to test Tool Shed installed
tools and only those tools.

::

    $ GALAXY_TEST_TOOL_CONF=config/tool_conf.xml sh run_tests.sh --report_file tool_tests_tool_conf.html functional.test_toolbox

The second command sets ``GALAXY_TEST_TOOL_CONF`` environment variable, which
will restrict the testing framework to considering a single tool conf file
(such as the default tools that ship with Galaxy
``config/tool_conf.xml.sample`` and which must have their dependencies setup
manually). The last argument to ``run_tests.sh``, ``functional.test_toolbox``
tells the test framework to run all the tool tests in the configured tool conf
file.

.. note:: *Tip:* To speed up tests you can use a pre-migrated database file the way Planemo
    does by setting the following environment variable before running
    ``run_tests.sh``.

    ::

         $ export GALAXY_TEST_DB_TEMPLATE="https://github.com/jmchilton/galaxy-downloads/raw/master/db_gx_rev_0127.sqlite"

.. _DOI: http://www.doi.org/
.. _BibTeX: http://www.bibtex.org/
.. _Dockerfile: https://docs.docker.com/reference/builder/
.. _Docker Hub: https://hub.docker.com/
.. _Planemo: http://planemo.readthedocs.org/
.. _Picard: http://broadinstitute.github.io/picard/

----------------------------------------------------------------------------
\.\.\. test tools against a package or container in a bioconda pull request?
----------------------------------------------------------------------------

First, obtain the artifacts of the PR by adding this comment:
``@BiocondaBot please fetch artifacts``. In the reply one finds the links to the
built package and docker image. 

In order to test the tool with the package add the following to the planemo call::

     $ planemo test ... --conda_channels LINK_TO_PACKAGE,conda-forge,bioconda,defaults ...

For containerized testing the docker image needs to be loaded::

     $ curl  -L "LINK_TO_DOCKER_IMAGE.tar.gz" | gzip -dc | docker load

A planemo test will then simply use this image::

     $ planemo test ... --biocontainers --no_conda_auto_init ...

--------------------------------------
\.\.\. interactively debug tool tests?
--------------------------------------

It can be desirable to interactively debug a tool test. In order to do so, start ``planemo test``
with the option ``--no_cleanup``. Inspect the output: After Galaxy starts up, the tests commence. At the
start of each test one finds a message: ``( <TOOL_ID> ) > Test-N``. After some upload jobs, the
actual tool job is started (it is the last before the next test is executed). There you will find
a message like ``Built script [/tmp/tmp1zixgse3/job_working_directory/000/3/tool_script.sh]``

In this case ``/tmp/tmp1zixgse3/job_working_directory/000/3/`` is the job dir. It contains some
files and directories of interest: 

- ``tool_script.sh``: the bash script generated from the tool's ``command`` and ``version_command``
  tags plus some boiler plate code
- ``galaxy_3.sh`` (note that the number may be different): a shell script setting up the environment
  (e.g. paths and environment variables), starting the ``tool_script.sh``, and postprocessing
  (e.g. error handling and setting metadata)
- ``working``: the job working directory
- ``outputs``: a directory containing the job stderr and stdout

For a tool test that uses a conda environment to resolve the requirements one can simply change
into ``working`` and execute ``../tool_script.sh`` (works as long as no special environment variables
are used; in this case ``../galaxy_3.sh`` needs to be executed after cleaning the job dir). 
By editing the tool script one may understand/fix problems in the ``command`` block faster than by
rerunning ``planemo test`` over and over again.

Alternatively one can change into the ``working`` dir and load the conda environment
(the code to do so can be found in ``tool_script.sh``: ``. PATH_TO_CONDA_ENV activate``). 
Afterwards one can execute individual commands, e.g. those found in ``tool_script.sh`` or variants.

For a tool test that uses Docker to to resolve the requirements one needs to execute 
``../galaxy_3.sh``, because it executes ``docker run ... tool_script.sh`` in order to rerun the job
(with a possible edited version of the tool script). In order to run the docker container 
interactively execute the ``docker run .... /bin/bash`` that you find in ``../galaxy_3.sh``
(i.e. ommitting the call of the ``tool_script.sh``) with added parameter ``-it``. 
