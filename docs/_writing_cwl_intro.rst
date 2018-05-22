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
flag tells Planemo to generate a Common Workflow Language tool. ``--id`` is
a short identifier for this tool and it should be unique across all tools.
``--name`` is a short, human-readable name for the the tool - it corresponds
to the ``label`` attribute in the CWL tool document.

The above command will generate the file ``seqtk_seq.cwl`` - which should look
like this.

.. literalinclude:: writing/seqtk_seq_v1.cwl
   :language: yaml

This tool file has the common fields required for a CWL tool with TODO notes,
but you will still need to open up the editor and fill out the command, describe
input parameters, tool outputs, writeup usage documentation (``doc``), etc..

The ``tool_init`` command can do a little bit better than this as well. We can
use the test command we tried above ``seqtk seq -A 2.fastq > 2.fasta`` as
an example to generate a command block by specifing the inputs and the outputs
as follows.

::

    $ planemo tool_init --force \
                        --cwl \
                        --id 'seqtk_seq' \
                        --name 'Convert to FASTA (seqtk)' \
                        --example_command 'seqtk seq -A 2.fastq > 2.fasta' \
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
                        --example_command 'seqtk seq -A 2.fastq > 2.fasta' \
                        --example_input 2.fastq \
                        --example_output 2.fasta \
                        --requirement seqtk@1.2 \
                        --container 'quay.io/biocontainers/seqtk:1.2--0' \
                        --test_case \
                        --help_from_command 'seqtk seq'

This command generates the following CWL YAML file.

.. literalinclude:: writing/seqtk_seq_v3.cwl
   :language: yaml

In addition to generating a CWL tool adding the ``--test_case`` flag generates from more files
that are useful including ``seqtk_seq_job.yml`` as shown below:

.. literalinclude:: writing/seqtk_seq_v3_job.yml
   :language: yaml

This is a CWL job input document and should allow you to run the example command using any CWL 
implementation. For instance if you have cwltool_ (``cwltool``) or Toil_ (``cwltoil``) on your
``PATH`` the following examples should work.

::

    $ cwltool seqtk_seq.cwl seqtk_seq_job.yml
    /Users/john/workspace/planemo/.venv/bin/cwltool 1.0.20180508202931
    Resolved 'seqtk_seq.cwl' to 'file:///Users/john/tool_init_exercise/seqtk_seq.cwl'
    [job seqtk_seq.cwl] /private/tmp/docker_tmpXgtSLt$ docker \
        run \
        -i \
        --volume=/private/tmp/docker_tmpXgtSLt:/private/var/spool/cwl:rw \
        --volume=/private/var/folders/78/zxz5mz4d0jn53xf0l06j7ppc0000gp/T/tmpGG1thW:/tmp:rw \
        --volume=/Users/john/tool_init_exercise/test-data/2.fastq:/private/var/lib/cwl/stg7db12d3a-2375-42ed-ba60-8a0ef69ffe80/2.fastq:ro \
        --workdir=/private/var/spool/cwl \
        --read-only=true \
        --log-driver=none \
        --user=502:20 \
        --rm \
        --env=TMPDIR=/tmp \
        --env=HOME=/private/var/spool/cwl \
        quay.io/biocontainers/seqtk:1.2--1 \
        seqtk \
        seq \
        -A \
        /private/var/lib/cwl/stg7db12d3a-2375-42ed-ba60-8a0ef69ffe80/2.fastq > /private/tmp/docker_tmpXgtSLt/out
    [job seqtk_seq.cwl] completed success
    {
        "output1": {
            "checksum": "sha1$322e001e5a99f19abdce9f02ad0f02a17b5066c2",
            "basename": "out",
            "location": "file:///Users/john/tool_init_exercise/out",
            "path": "/Users/john/tool_init_exercise/out",
            "class": "File",
            "size": 150
        }
    }

::

    $ cwltoil seqtk_seq.cwl seqtk_seq_job.yml
    jlaptop17.local 2018-05-21 15:25:30,630 MainThread INFO toil.lib.bioio: Root logger is at level 'INFO', 'toil' logger at level 'INFO'.
    jlaptop17.local 2018-05-21 15:25:30,648 MainThread INFO toil.jobStores.abstractJobStore: The workflow ID is: '55a08d91-1852-4069-97a9-741abd2ea04e'
    Resolved 'seqtk_seq.cwl' to 'file:///Users/john/tool_init_exercise/seqtk_seq.cwl'
    jlaptop17.local 2018-05-21 15:25:30,650 MainThread INFO cwltool: Resolved 'seqtk_seq.cwl' to 'file:///Users/john/tool_init_exercise/seqtk_seq.cwl'
    jlaptop17.local 2018-05-21 15:25:31,793 MainThread INFO toil.common: Using the single machine batch system
    jlaptop17.local 2018-05-21 15:25:31,793 MainThread WARNING toil.batchSystems.singleMachine: Limiting maxCores to CPU count of system (8).
    jlaptop17.local 2018-05-21 15:25:31,793 MainThread WARNING toil.batchSystems.singleMachine: Limiting maxMemory to physically available memory (17179869184).
    jlaptop17.local 2018-05-21 15:25:31,800 MainThread INFO toil.common: Created the workflow directory at /var/folders/78/zxz5mz4d0jn53xf0l06j7ppc0000gp/T/toil-55a08d91-1852-4069-97a9-741abd2ea04e-132281828025877
    jlaptop17.local 2018-05-21 15:25:31,800 MainThread WARNING toil.batchSystems.singleMachine: Limiting maxDisk to physically available disk (206962089984).
    jlaptop17.local 2018-05-21 15:25:31,808 MainThread INFO toil.common: User script ModuleDescriptor(dirPath='/Users/john/workspace/planemo/.venv/lib/python2.7/site-packages', name='toil.cwl.cwltoil', fromVirtualEnv=True) belongs to Toil. No need to auto-deploy it.
    jlaptop17.local 2018-05-21 15:25:31,809 MainThread INFO toil.common: No user script to auto-deploy.
    jlaptop17.local 2018-05-21 15:25:31,809 MainThread INFO toil.common: Written the environment for the jobs to the environment file
    jlaptop17.local 2018-05-21 15:25:31,809 MainThread INFO toil.common: Caching all jobs in job store
    jlaptop17.local 2018-05-21 15:25:31,809 MainThread INFO toil.common: 0 jobs downloaded.
    jlaptop17.local 2018-05-21 15:25:31,825 MainThread INFO toil: Running Toil version 3.15.0-0e3a87e738f5e0e7cff64bfdad337d592bd92704.
    jlaptop17.local 2018-05-21 15:25:31,825 MainThread INFO toil.realtimeLogger: Real-time logging disabled
    jlaptop17.local 2018-05-21 15:25:31,832 MainThread INFO toil.toilState: (Re)building internal scheduler state
    2018-05-21 15:25:31,832 - toil.toilState - INFO - (Re)building internal scheduler state
    jlaptop17.local 2018-05-21 15:25:31,832 MainThread INFO toil.leader: Found 1 jobs to start and 0 jobs with successors to run
    2018-05-21 15:25:31,832 - toil.leader - INFO - Found 1 jobs to start and 0 jobs with successors to run
    jlaptop17.local 2018-05-21 15:25:31,832 MainThread INFO toil.leader: Checked batch system has no running jobs and no updated jobs
    2018-05-21 15:25:31,832 - toil.leader - INFO - Checked batch system has no running jobs and no updated jobs
    jlaptop17.local 2018-05-21 15:25:31,833 MainThread INFO toil.leader: Starting the main loop
    2018-05-21 15:25:31,833 - toil.leader - INFO - Starting the main loop
    jlaptop17.local 2018-05-21 15:25:31,834 MainThread INFO toil.leader: Issued job 'file:///Users/john/tool_init_exercise/seqtk_seq.cwl' seqtk seq u/c/jobzYKQ3V with job batch system ID: 0 and cores: 1, disk: 3.0 G, and memory: 2.0 G
    2018-05-21 15:25:31,834 - toil.leader - INFO - Issued job 'file:///Users/john/tool_init_exercise/seqtk_seq.cwl' seqtk seq u/c/jobzYKQ3V with job batch system ID: 0 and cores: 1, disk: 3.0 G, and memory: 2.0 G
    jlaptop17.local 2018-05-21 15:25:33,953 MainThread INFO toil.leader: Job ended successfully: 'file:///Users/john/tool_init_exercise/seqtk_seq.cwl' seqtk seq u/c/jobzYKQ3V
    2018-05-21 15:25:33,953 - toil.leader - INFO - Job ended successfully: 'file:///Users/john/tool_init_exercise/seqtk_seq.cwl' seqtk seq u/c/jobzYKQ3V
    jlaptop17.local 2018-05-21 15:25:33,955 MainThread INFO toil.leader: Finished the main loop: no jobs left to run
    2018-05-21 15:25:33,955 - toil.leader - INFO - Finished the main loop: no jobs left to run
    jlaptop17.local 2018-05-21 15:25:33,955 MainThread INFO toil.serviceManager: Waiting for service manager thread to finish ...
    2018-05-21 15:25:33,955 - toil.serviceManager - INFO - Waiting for service manager thread to finish ...
    jlaptop17.local 2018-05-21 15:25:34,841 MainThread INFO toil.serviceManager: ... finished shutting down the service manager. Took 0.885795116425 seconds
    2018-05-21 15:25:34,841 - toil.serviceManager - INFO - ... finished shutting down the service manager. Took 0.885795116425 seconds
    jlaptop17.local 2018-05-21 15:25:34,842 MainThread INFO toil.statsAndLogging: Waiting for stats and logging collator thread to finish ...
    2018-05-21 15:25:34,842 - toil.statsAndLogging - INFO - Waiting for stats and logging collator thread to finish ...
    jlaptop17.local 2018-05-21 15:25:34,854 MainThread INFO toil.statsAndLogging: ... finished collating stats and logs. Took 0.0120511054993 seconds
    2018-05-21 15:25:34,854 - toil.statsAndLogging - INFO - ... finished collating stats and logs. Took 0.0120511054993 seconds
    jlaptop17.local 2018-05-21 15:25:34,855 MainThread INFO toil.leader: Finished toil run successfully
    2018-05-21 15:25:34,855 - toil.leader - INFO - Finished toil run successfully
    {
        "output1": {
            "checksum": "sha1$322e001e5a99f19abdce9f02ad0f02a17b5066c2",
            "basename": "out",
            "nameext": "",
            "nameroot": "out",
            "http://commonwl.org/cwltool#generation": 0,
            "location": "file:///Users/john/tool_init_exercise/out",
            "class": "File",
            "size": 150
        }
    jlaptop17.local 2018-05-21 15:25:34,866 MainThread INFO toil.common: Successfully deleted the job store: <toil.jobStores.fileJobStore.FileJobStore object at 0x1057205d0>
    }2018-05-21 15:25:34,866 - toil.common - INFO - Successfully deleted the job store: <toil.jobStores.fileJobStore.FileJobStore object at 0x1057205d0>


At this point we have a fairly a functional CWL tool with test and usage
documentation. This was a pretty simple example - usually you will need to
put more work into the tool to get to this point - ``tool_init`` is really
just designed to get you started.

Now lets lint and test the tool we have developed. The Planemo's ``lint`` (or
just ``l``) command will review tool for validity, obvious mistakes, and
Planemo "best practices".

::

    $ planemo l seqtk_seq.cwl
    Linting tool /Users/john/workspace/planemo/docs/writing/seqtk_seq.cwl
    Applying linter general... CHECK
    .. CHECK: Tool defines a version [0.0.1].
    .. CHECK: Tool defines a name [Convert to FASTA (seqtk)].
    .. CHECK: Tool defines an id [seqtk_seq_v3].
    .. CHECK: Tool specifies profile version [16.04].
    Applying linter cwl_validation... CHECK
    .. INFO: CWL appears to be valid.
    Applying linter docker_image... CHECK
    .. INFO: Tool will run in Docker image [quay.io/biocontainers/seqtk:1.2--1].
    Applying linter new_draft... CHECK
    .. INFO: Modern CWL version [v1.0]

In addition to the actual tool and job files, ``--test_case`` caused a test file
to be generated using the example command and provided test data. The file contents
are as follows:

.. literalinclude:: writing/seqtk_seq_v3_tests.yml
   :language: yaml

Unlike the job file, this file is a Planemo-specific artifact. This file may contain
1 or more tests - each test is an element of the top-level list. ``tool_init`` will use
the example command to build just one test.

Each test consists of a few parts:

- ``doc`` - this attribute provides a short description for the test.
- ``job`` - this can be the path to a CWL job description or a job
  description embedded right in the test (``tool_init`` builds the latter).
- ``outputs`` - this section describes the expected output for a test. Each
  output ID of the tool or workflow under test can appear as a key. The
  example above just describes expected specific output file contents exactly
  but many more expectations can be described.

For more information on the test file format check out `the Test Format docs
<http://planemo.readthedocs.io/en/latest/test_format.html>`__.

The tests described in this file can be run using the ``planemo test`` command
on the original file.

::

    $ planemo test --no-container seqtk_seq.cwl
    Enable beta testing mode for testing.
    cwltool INFO: /Users/john/workspace/planemo/.venv/bin/planemo 1.0.20180508202931
    cwltool INFO: Resolved '/Users/john/tool_init_exercise/seqtk_seq.cwl' to 'file:///Users/john/tool_init_exercise/seqtk_seq.cwl'
    cwltool INFO: [job seqtk_seq.cwl] /private/tmp/docker_tmpvLE9SS$ seqtk \
        seq \
        -A \
        /private/var/folders/78/zxz5mz4d0jn53xf0l06j7ppc0000gp/T/tmpGM22d_/stg0c0cad75-7ca0-4f3a-9d77-63e9c49f5353/2.fastq > /private/tmp/docker_tmpvLE9SS/out
    cwltool INFO: [job seqtk_seq.cwl] completed success
    cwltool INFO: Final process status is success
    All 1 test(s) executed passed.
    seqtk_seq_0: passed

This is a bit a different than running the job. For one thing, we don't need to
specify an input job - instead Planemo will automatically find the test file and
run all the jobs described inside that file. Additionally, Planemo will check the
outputs to ensure the match the test expectations.

.. include:: _writing_test_reports.rst

The above test example used cwltool_ to run our test and disabled containerization.
By dropping the ``--no-container`` argument we can run the tool in a Docker container.
By passing an engine argument as ``--engine toil`` we can run our test in Toil_, an
alternative CWL implementation.

::

    $ planemo test seqtk_seq.cwl
    Enable beta testing mode for testing.
    cwltool INFO: /Users/john/workspace/planemo/.venv/bin/planemo 1.0.20180508202931
    cwltool INFO: Resolved '/Users/john/tool_init_exercise/seqtk_seq.cwl' to 'file:///Users/john/tool_init_exercise/seqtk_seq.cwl'
    cwltool INFO: [job seqtk_seq.cwl] /private/tmp/docker_tmpUeIpXJ$ docker \
        run \
        -i \
        --volume=/private/tmp/docker_tmpUeIpXJ:/private/var/spool/cwl:rw \
        --volume=/private/var/folders/78/zxz5mz4d0jn53xf0l06j7ppc0000gp/T/tmpteo_2Z:/tmp:rw \
        --volume=/Users/john/tool_init_exercise/test-data/2.fastq:/private/var/lib/cwl/stg939ee60b-a194-4177-8410-c40a1acb38ea/2.fastq:ro \
        --workdir=/private/var/spool/cwl \
        --read-only=true \
        --log-driver=none \
        --user=502:20 \
        --rm \
        --env=TMPDIR=/tmp \
        --env=HOME=/private/var/spool/cwl \
        quay.io/biocontainers/seqtk:1.2--1 \
        seqtk \
        seq \
        -A \
        /private/var/lib/cwl/stg939ee60b-a194-4177-8410-c40a1acb38ea/2.fastq > /private/tmp/docker_tmpUeIpXJ/out
    cwltool INFO: [job seqtk_seq.cwl] completed success
    cwltool INFO: Final process status is success
    All 1 test(s) executed passed.
    seqtk_seq_0: passed
    $ planemo test --no-container --engine toil seqtk_seq.cwl
    Enable beta testing mode for testing.
    All 1 test(s) executed passed.
    seqtk_seq_0: passed

For more information on the Common Workflow Language check out the Draft 3
`User Guide`_ and Specification_.

.. _YAML: http://yaml.org/
.. _User Guide: http://www.commonwl.org/user_guide/
.. _Specification: http://www.commonwl.org/v1.0/CommandLineTool.html
.. _Toil: https://github.com/BD2KGenomics/toil
.. _cwltool: https://github.com/common-workflow-language/cwltool
