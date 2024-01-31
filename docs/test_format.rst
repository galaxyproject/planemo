====================
Test Format
====================

Planemo has traditionally been used to test Galaxy_ tools.

::

    $ planemo test galaxy_tool.xml

This starts a Galaxy instance, runs the tests described in the XML file, prints a nice summary of the
test results (pass or fail for each test) in the console and creates an HTML report in the current
directory. Additional bells and whistles include the ability to generate XUnit reports, publish
test results and get embedded Markdown to link to them for PRs, and test remote artifacts in Git repositories.

For more information about testing Galaxy tools using embedded tool XML tests see the tutorial-style chapter
`Test-Driven Development <https://planemo.readthedocs.io/en/latest/writing_advanced.html#test-driven-development>`__
of Galaxy tools.

Much of this same functionality is now also available for Galaxy_ Workflows as well as `Common Workflow Language`_
(CWL) tools and workflows. The rest of this page describes the test format and testing options for these
artifacts.

Unlike the traditional Galaxy tool approach, these newer types of artifacts should define tests in files
located next to the artifact. For instance, if ``planemo test`` is called on a Galaxy workflow called ``ref-rnaseq.ga``
tests should be defined in ``ref-rnaseq-tests.yml`` or ``ref-rnaseq-tests.yaml``. If instead it is called on a
CWL_ tool called ``seqtk_seq.cwl``, tests can be defined in ``seqtk_seq_tests.yml`` for instance.

Below are two examples of such YAML files - the first for a CWL_ tool and the second for Galaxy_ workflow. Note the
same testing file format is used for both kinds of artifacts.

.. code-block:: yaml

    - doc: simple usage test
      job: pear_job.yml
      outputs:
        assembled_pairs:
          path: test-data/pear_assembled_results1.fastq
        unassembled_forward_reads:
          path: test-data/pear_unassembled_forward_results1.fastq


.. code-block:: yaml

    - doc: Test sample data for Microbial variant calling workflow
      job:
        mutant_R1:
          class: File
          path: mutant_R1.fastq
        mutant_R2:
          class: File
          path: mutant_R2.fastq
        wildtype.fna:
          class: File
          location: https://zenodo.org/record/582600/files/wildtype.fna
        wildtype.gbk:
          class: File
          location: https://zenodo.org/record/582600/files/wildtype.gbk
        wildtype.gff:
          class: File
          location: https://zenodo.org/record/582600/files/wildtype.gff
      outputs:
        jbrowse_html:
          asserts:
            has_text:
              text: "JBrowseDefaultMainPage"
        snippy_fasta:
          asserts:
            has_line:
              line: '>Wildtype Staphylococcus aureus strain WT.'
        snippy_tabular:
          asserts:
            has_n_columns:
              n: 2

The above examples illustrate that each test file is broken into a list of test
cases. Each test case should have a ``doc`` describing the test, a ``job`` description
the describes the inputs for an execution of the target artifact, and an ``outputs``
mapping that describes assertions about outputs to test.


``job``
--------

The ``job`` object can be a mapping embedded right in the test file or a reference to a
an external "job" input file. The job input file is a proper CWL_ job document - which
is fairly straight forward as demonstrated in the above examples. Planemo adapts the
CWL_ job document to Galaxy_ workflows and tools - using input names for Galaxy_ tools
and input node labels for workflows.

Input files can be specified using either ``path`` attributes (which should generally be file
paths relative to the artifact and test directory) or ``location`` (which should be a URI). The
examples above demonstrate using both paths relative to the tool file and test data published
to `Zenodo <https://zenodo.org/>`__.

Embedded job objects result in cleaner test suites that are simpler to read. One advantage of
instead using external job input files is that the job object can be reused to invoke the
runnable artifact outside the context of testing with ``planemo run``.

.. note::

    These job objects can be run directly with ``planemo run``.

    ::

        $ planemo run --engine=<engine_type> [ENGINE_OPTIONS] [ARTIFACT_PATH] [JOB_PATH]

    This should be familar to CWL developers - and indeed with ``--engine=cwltool`` this works as a formal CWL
    runner. Planemo provides a uniform interface to Galaxy for Galaxy workflows and tools though using the same
    CLI invocation if ``--engine=galaxy`` (for a Planemo managed Galaxy instance), ``--engine=docker_galaxy``
    (for a Docker instance of Galaxy launched by Planemo), or ``--engine=external_galaxy`` (for a running
    remote Galaxy instance).

Certain Galaxy objects don't map cleanly to CWL_ job objects so Planemo attempts to extend
the format with new constructs for running and testing Galaxy objects - such as describing 
collections and composite inputs.

Galaxy Collection Inputs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following example demonstrates two ways to create input lists for Galaxy tests.

.. literalinclude:: test_example_collection_input.yml
   :language: yaml

Simply specifying files in YAML lists in the input job (like vanilla CWL_
job descriptions) will result in a simple Galaxy list. This is simple but the downside is you
have no control of the list identifiers - which are often important in Galaxy workflows. When
more control is desired, you may describe an explicit Galaxy collection with an input object of
``class: Collection``. This variant (also shown in the above example) allows creating collections
of type other than ``list`` and allows specifying element identifiers with the ``identifier``
declaration under the list of collection ``elements``.

The explicit Galaxy collection creation syntax also makes describing nested collections such as
lists of pairs very natural. The following example is used in Planemo's test suite to illustrate
this:

.. literalinclude:: test_example_nested_collection_input.yml
   :language: yaml

Galaxy Composite Inputs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The syntax for specifying composite inputs is a little more basic still and simply must be specified
as a list of local files (mirroring Galaxy Tool XML test syntax). While ``class`` is assumed to
be ``File`` and URIs aren't yet tested.


.. literalinclude:: test_example_composite_input.yml
   :language: yaml

Galaxy Tags
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Requires Galaxy 20.09 or newer.**

Tags and group tags play important roles in many Galaxy workflows. These can be tested
by simply add a list of ``tags:`` to the YAML corresponding to the dataset in the
collection. The following example demonstrates this:

.. literalinclude:: test_example_tagged_input.yml
   :language: yaml

``outputs``
--------------

Galaxy_ tools and CWL_ artifacts have obvious output names that much match the mapping in this block on test
file. Galaxy workflows require explicit output labels to be used with tests, but the important outputs in
your workflows should be labeled anyway to work with Galaxy subworkflows and more cleanly with API calls.

If an output is known, fixed, and small it makes a lot of sense to just include a copy of the output next
to your test and set ``file: relative/path/to/output`` in your output definition block as show in the first
example above. For completely reproducible processes this is a great guarantee that results are fixed over
time, across CWL_ engines and engine versions. If the results are fixed but large - it may make sense to just
describe the outputs by a SHA1_ checksum_.

.. code-block:: yaml

    - doc: Simple concat workflow test
      job: wf1.gxwf-job.yml
      outputs:
        wf_output_1:
          checksum: "sha1$a0b65939670bc2c010f4d5d6a0b3e4e4590fb92b"

One advantage of including an exact file instead of a checksum is that Planemo can produce very nice line
by line diffs for incorrect test results by comparing an expected output to an actual output.

There are reasons one may not be able to write exact test assertions about outputs however.
Perhaps date or time information is incorporated into a result, unseeded random numbers are used, small numeric
differences occur across runtimes of interest, etc..
Even if one can write exact tests, a really useful technique is to write more liberal sanity checks on outputs as one
builds up workflows that may be changing rapidly and develops complex tools or workflows via a
`Test-Driven Development cycle
<https://en.wikipedia.org/wiki/Test-driven_development#Test-driven_development_cycle>`__
using Planemo. *Tests shouldn't just be an extra step you have to do after development is done, they should
guide development as well.*

In all of these cases, a variety of other assertions can be run against the execution results to verify outputs.
The "Microbial variant calling workflow" example at the beginning of this chapter demonstrates some assertions one can make about the contents of result files.
Some additional examples of inexact file comparisons taken from an artificial test case in the Planemo test suite are shown below.

.. literalinclude:: example_assertions.yml
   :language: yaml

Currently, the full list of available assertions is only documented as part of the `Galaxy Tool XML format <https://docs.galaxyproject.org/en/latest/dev/schema.html>`__ definition in the section on `asserting the contents of Galaxy tool outputs <https://docs.galaxyproject.org/en/latest/dev/schema.html#tool-tests-test-output-assert-contents>`__, but it should be fairly easy to translate this XML syntax into the YAML format above.


Engines for Testing
---------------------

Below are descriptions of various testing engines that can be used with Planemo (both with the
`test command`_ and the `run command`_) as well as some command-line options of particular interest for testing.
The first two types ``cwltool`` and ``toil`` can be used to test CWL artifacts (tools and workflows).
The remaining engine types are variations on engines that target Galaxy and are useful for testing
workflows (and tools with newer style tests or job documents).

``cwltool``
~~~~~~~~~~~~~~~~~

::

    $ planemo test --engine cwltool [--no-container] [--biocontainers]

This is the most straight forward engine, it can be used to test CWL_ tools and workflows using
the CWL reference implementation cwltool_ (bundled as a dependency of Planemo). Use the ``--no-container``
option to disable Docker and use Conda resolution of ``SoftwareRequirement``s or applications
on the ``PATH``. Use the ``--biocontainers`` flag to use BioContainers_ for tools without
explicit ``DockerRequirement`` hints.

``toil``
~~~~~~~~~~~~~~~~~

::

    $ planemo test --engine toil [--no-container] [--biocontainers]

This engine largely mirrors the ``cwltool`` engine but runs CWL artifacts using Toil_. Toil_ is
an optional dependency of Planemo so you will likely have to install it in Planemo's environment
using ``pip install toil``.

``galaxy``
~~~~~~~~~~~~~~~~~

::

    $ planemo test [--docker] [--biocontainers] [--profile <profile>] [--galaxy_root <path>] [--extra_tools <path>]

This is the default engine type, but can be made explicit ``--engine galaxy``. With this engine
Planemo will start a Galaxy instance and test against it.

Planemo will automatically detect and load "stock" Galaxy tools used by workflows and install any
Tool Shed tools contained in the workflow, if other non-Tool Shed tools are required for a workflow they can be
loaded using ``--extra_tools``.

Set ``--galaxy_root`` to target an externally cloned Galaxy directory or use ``--galaxy_branch`` to
target a particular branch of Galaxy other than the latest stable.

Use the ``--biocontainers`` flag to enable Docker and use BioContainers_ for tools or use ``--docker`` to use
Docker but limited to tools configured with ``container`` tags.

By default Galaxy when configured by Planemo will attempt to run with an sqlite database. This configuration
is quite buggy and should not be used to test workflows. The ``--profile`` option can be used to target
a pre-configured Postgres database created with ``planemo profile_create`` and use it for testing. In
addition to making Galaxy more robust this should speed up testing after the initial setup of the database.

::

    planemo profile_create --database_type [postgres|docker_postgres] my_cool_name
    planemo test --profile my_cool_name

If ``--database_type`` is specified as ``docker_postgres``, Planemo will attempt to startup a postgres server
in a Docker container automatically for testing. If instead ``postgres`` is specified Planemo will attempt
to interact with Postgres using ``psql`` (assumed to be on the ``PATH``). For a description on more Postgres
connection options check out the documentation for the ``database_create`` `command
<http://planemo.readthedocs.io/en/latest/commands.html#database-create-command>`__ that has similar options.

Profiles may also really help testing local setups by saving previously installed shed repository installations
and Conda environments.

``docker_galaxy``
~~~~~~~~~~~~~~~~~

::

    $ planemo test --engine docker_galaxy [--extra_tools <path>] [--docker_extra_volume <path>] [--docker_galaxy_image <image>]

With this engine Planemo will start a Docker container to run tests against it. See the `docker-galaxy-stable
<https://github.com/bgruening/docker-galaxy-stable/>`__ project spearheaded by Björn Grüning for more
information on Docker-ized Galaxy execution. The exact container image to use can be controlled using the
``--docker_galaxy_image`` option.

Planemo will automatically detect and load "stock" Galaxy tools used by workflows and install any
Tool Shed tools contained in the workflow, if other non-Tool Shed tools are required - they can be
loaded using ``--extra_tools``.

At the time of this writing, there is a bug in Planemo that requires using the ``--docker_extra_volume``
option to mount test data into the testing container.

``external_galaxy``
~~~~~~~~~~~~~~~~~~~~

    $ planemo test --engine external_galaxy --galaxy_admin_key <admin_key> --galaxy_user_key <user_key> [--no_shed_install] [--polling_backoff <integer>] --galaxy_url <url>

This is primarily useful for testing workflows against already running Galaxy instances. An admin or
bootstrap API key should be supplied to install missing tool repositories for the workflow and a user API
key should be supplied to run the workflow using. If you wish to skip tool shed repository installation
(this requires all the tools be present already), use the ``--no_shed_install`` option. If you want to
reduce the load on the target Galaxy while checking for the status changes use the ``--polling_backoff <integer>``
option where integer is the incremental increase in seconds for every request.

To run tool tests against a running Galaxy, ``galaxy-tool-test`` is a script that gets installed with
galaxy-tool-util and so may very well already be on your ``PATH``. Check out the options available with that
using ``galaxy-tool-test --help``. If you're interested in running all the tool tests corresponding to
a workflow on a running server, check out the galaxy-workflow-tool-tests_ project that is a wrapper around
galaxy-tool-test that has all the same options but that filters to the tool tests to just run those
from a specific workflow.

This engine can also be used to test workflows already available in the running Galaxy instance.
While you don't need to download and synchronize the target workflow on your local filesystem, you
do need to provide a path to find the test definition and test data paths.

An example of doing this is included in Planemo's `test data <https://github.com/galaxyproject/planemo/tree/master/tests/data>`__.
The workflow test definition ``wf11-remote.gxwf-test.yml`` exists but no corresponding workflow file ``wf11-remote.gxwf.yml`` exists.
The workflow is assumed to already exist in some Galaxy server. For instance, it might exist somewhere with id ``99113b2b119318e1``.
Then ``planemo test`` could be run with ``gxid://workflows/99113b2b119318e1?runnable_path=/path/to/wf11-remote.gxwf.yml``
as the last argument to test this workflow with that test data. Note this path ``/path/to/wf11-remote.gxwf.yml``
doesn't need to exist, but it is used to find ``wf11-remote.gxwf-test.yml``.

Galaxy Testing Template
-------------------------

The following is a script that can be used with `continuous integration`_ (CI) services such
Travis_ to test Galaxy workflows in a Github repository. This shell script can be configured via
various environment variables and shows off some of the modalities Planemo ``test`` should work in
(there may be bugs but we are trying to stablize this functionality).

.. literalinclude:: run_galaxy_workflow_tests.sh
   :language: sh

A Travis_ configuration file (.travis.yml) that would test workflows using a Docker Galaxy image might
look like:

.. literalinclude:: example_travis_gxwf_docker.yml
   :language: yaml

To skip Docker and instead test with a native Galaxy instance and postgres database one might use the
configuration:

.. literalinclude:: example_travis_gxwf_native.yml
   :language: yaml

.. _Travis: https://travis-ci.org/
.. _test command: http://planemo.readthedocs.org/en/latest/commands.html#test-command
.. _run command: http://planemo.readthedocs.org/en/latest/commands.html#run-command
.. _SHA1: https://en.wikipedia.org/wiki/SHA-1
.. _checksum: https://en.wikipedia.org/wiki/Checksum
.. _continuous integration: https://en.wikipedia.org/wiki/Continuous_integration
.. _Common Workflow Language: https://www.commonwl.org/
.. _CWL: https://www.commonwl.org/
.. _cwltool: https://github.com/common-workflow-language/cwltool
.. _Galaxy: http://galaxyproject.org/
.. _Toil: https://github.com/BD2KGenomics/toil
.. _BioContainers: http://biocontainers.pro/
.. _galaxy-workflow-tool-tests: https://github.com/galaxy-works/galaxy-workflow-tool-tests
