.. _dependencies_and_conda_cwl:

Dependencies and Conda
===========================================

----------------------------------------------------------------
Specifying and Using `Software Requirements`_
----------------------------------------------------------------

.. include:: _writing_conda_init.rst

.. note:: Why not just use containers?

    Containers are great, use containers (be it Docker_, Singularity_, etc.) whenever possible to
    increase reproducibility and portability of your tools and workflow. Building ad hoc containers
    to support CWL tools (e.g. custom ``Dockerfile`` definitions) has serious limitations, in the next
    tutorial on containers we will argue that using Biocontainers_ built or discovered
    from your tool's `Software Requirements`_ is a superior approach.

    Besides leading to better containers, there are other reasons to describe
    `Software Requirements`_ also - it will allow your tool to be used in environments without
    container runtimes available and provides valuable and actionable metadata about the computation
    described by the tool.

    Read more about this whole dependency stack in our preprint `Practical computational reproducibility
    in the life sciences <https://www.biorxiv.org/content/early/2017/10/10/200683>`__

The `Common Workflow Language`_ specification loosely describes
`Software Requirements`_ - a way to map CWL hints to packages, environment
modules, or any other mechanism to describe dependencies for running a tool
outside of a container. The large and active Galaxy tool development community
has built an open source library and set of best practices for describing dependencies
for Galaxy that should work just as well for CWL. The library has been integrated
with cwltool_ and Toil_ to enable CWL tool authors and users to leverage the
power and flexibility of the Galaxy dependency management and best practices.

While `Software Requirements`_ can be configured to resolve dependencies various ways,
Planemo is configured with opinionated defaults geared at making building CWL tools that
target Conda_ as easy as possible and build tools with requirements compatible with
cwltool_ and Toil_ when running outside containers.

During the tool development `introductory tutorial`_, we called ``planemo tool_init``
with the argument ``--requirement seqtk@1.2`` and the resulting tool contained such
a ``SoftwareRequirement`` in the form the following the YAML fragment::

    SoftwareRequirement:
      packages:
      - package: seqtk
        version:
        - "1.2"

Planemo (and cwltool_ and Toil_) can interpret these ``SoftwareRequirement`` annotations in various ways
including as Conda packages. When interpreting these as Conda packages
these runtimes can setup isolated, reproducible Conda environments for tool execution with the correct
packages installed (e.g. ``seqtk`` in the above example).

.. include:: _writing_conda_overview.rst

We can check if the requirements on a tool are available in best practice
Conda channels using an extended form of the ``planemo lint`` command (``planemo lint`` was
introduced in the `introductory tutorial`_). Passing ``--conda_requirements`` flag will ensure all
listed requirements are found.

::

    $ planemo lint --conda_requirements seqtk_seq.cwl
    Linting tool /Users/john/workspace/planemo/docs/writing/seqtk_seq.cwl
      ...
    Applying linter requirements_in_conda... CHECK
    .. INFO: Requirement [seqtk@1.2] matches target in best practice Conda channel [https://conda.anaconda.org/bioconda/osx-64].

.. note:: You can download a more complete version of the CWL seqtk_ seq from the Planemo
    tutorial using the command::

        $ planemo project_init --template=seqtk_complete_cwl seqtk_example
        $ cd seqtk_example

We can verify these tool requirements install with the ``conda_install`` command. With
its default parameters ``conda_install`` processes tools and creates isolated environments
for their declared `Software Requirements`_ (mirroring what can be done in production with 
cwltool_ and Toil_).

::

    $ planemo conda_install seqtk_seq.cwl
    Install conda target CondaTarget[seqtk,version=1.2]
    /home/john/miniconda3/bin/conda create -y --name __seqtk@1.2 seqtk=1.2
    Fetching package metadata ...............
    Solving package specifications: ..........

    Package plan for installation in environment /home/john/miniconda2/envs/__seqtk@1.2:

    The following packages will be downloaded:

        package                    |            build
        ---------------------------|-----------------
        seqtk-1.2                  |                0          29 KB  bioconda

    The following NEW packages will be INSTALLED:

        seqtk: 1.2-0   bioconda
        zlib:  1.2.8-3

    Fetching packages ...
    seqtk-1.2-0.ta 100% |#############################################################| Time: 0:00:00 444.71 kB/s
    Extracting packages ...
    [      COMPLETE      ]|################################################################################| 100%
    Linking packages ...
    [      COMPLETE      ]|################################################################################| 100%
    #
    # To activate this environment, use:
    # > source activate __seqtk@1.2
    #
    # To deactivate this environment, use:
    # > source deactivate __seqtk@1.2
    #
    $ which seqtk
    seqtk not found
    $

The above install worked properly, but ``seqtk`` is not on your ``PATH`` because this merely
created an environment within the Conda directory for the seqtk installation. Planemo
will configure cwltool during testing to reuse this environment. If you wish to interactively explore
the resulting enviornment to explore the installed tool or produce test data the output
of the ``conda_env`` command can be sourced.

::

    $ . <(planemo conda_env seqtk_seq.cwl)
    Deactivate environment with conda_env_deactivate
    (seqtk_seq) $ which seqtk
    /home/planemo/miniconda3/envs/jobdepsiJClEUfecc6d406196737781ff4456ec60975c137e04884e4f4b05dc68192f7cec4656/bin/seqtk
    (seqtk_seq) $ seqtk seq

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
             -U        convert all bases to uppercases
             -S        strip of white spaces in sequences
    (seqtk_seq) $ conda_env_deactivate
    $

As shown above the ``conda_env_deactivate`` will be created in this environment and can
be used to restore your initial shell configuration.

Here is a portion of the output from the testing command ``planemo test seqtk_seq.cwl``
demonstrating using this tool.

::

    $ planemo test --no-container seqtk_seq.cwl
    Enable beta testing mode for testing.
    cwltool INFO: /Users/john/workspace/planemo/.venv/bin/planemo 1.0.20170828135420
    cwltool INFO: Resolved '/Users/john/workspace/planemo/project_templates/seqtk_complete_cwl/seqtk_seq.cwl' to 'file:///Users/john/workspace/planemo/project_templates/seqtk_complete_cwl/seqtk_seq.cwl'
    cwltool INFO: [job seqtk_seq.cwl] /private/var/folders/78/zxz5mz4d0jn53xf0l06j7ppc0000gp/T/tmpaDQ1nK$ seqtk \
        seq \
        -a \
        /private/var/folders/78/zxz5mz4d0jn53xf0l06j7ppc0000gp/T/tmpJtPKCr/stg24cf7e67-5ca6-44a4-a46b-26cbe104e1d4/2.fastq > /private/var/folders/78/zxz5mz4d0jn53xf0l06j7ppc0000gp/T/tmpaDQ1nK/out
    cwltool INFO: [job seqtk_seq.cwl] completed success
    cwltool INFO: Final process status is success
    galaxy.tools.parser.factory INFO: Loading CWL tool - this is experimental - tool likely will not function in future at least in same way.
    All 1 test(s) executed passed.
    seqtk_seq_0: passed

Since ``seqtk`` isn't on the path and we did not use a container, we can see the SoftwareRequirement 
resolution was successful and it found the environment we previously installed with ``conda_install``.

This can be used outside of Planemo testing as well, the following invocation shows running a job
with cwltool_ using an environment like the one created above:

::

    $ cwltool --no-container --beta-conda-dependencies seqtk_seq.cwl seqtk_seq_job.yml
    /Users/john/workspace/planemo/.venv/bin/cwltool 1.0.20180508202931
    Resolved 'seqtk_seq.cwl' to 'file:///Users/john/workspace/planemo/project_templates/seqtk_complete_cwl/seqtk_seq.cwl'
    No handlers could be found for logger "rdflib.term"
    [job seqtk_seq.cwl] /private/tmp/docker_tmpDQYeqC$ seqtk \
        seq \
        -a \
        /private/var/folders/78/zxz5mz4d0jn53xf0l06j7ppc0000gp/T/tmpQwBqPo/stg8cf2282a-d807-4f90-b94d-feeda004cacd/2.fastq > /private/tmp/docker_tmpDQYeqC/out
    PREFIX=/Users/john/workspace/planemo/project_templates/seqtk_complete_cwl/cwltool_deps/_conda
    installing: python-3.6.3-h47c878a_7 ...
    Python 3.6.3 :: Anaconda, Inc.
    installing: ca-certificates-2017.08.26-ha1e5d58_0 ...
    installing: conda-env-2.6.0-h36134e3_0 ...
    installing: libcxxabi-4.0.1-hebd6815_0 ...
    installing: tk-8.6.7-h35a86e2_3 ...
    installing: xz-5.2.3-h0278029_2 ...
    installing: yaml-0.1.7-hc338f04_2 ...
    installing: zlib-1.2.11-hf3cbc9b_2 ...
    installing: libcxx-4.0.1-h579ed51_0 ...
    installing: openssl-1.0.2n-hdbc3d79_0 ...
    installing: libffi-3.2.1-h475c297_4 ...
    installing: ncurses-6.0-hd04f020_2 ...
    installing: libedit-3.1-hb4e282d_0 ...
    installing: readline-7.0-hc1231fa_4 ...
    installing: sqlite-3.20.1-h7e4c145_2 ...
    installing: asn1crypto-0.23.0-py36h782d450_0 ...
    installing: certifi-2017.11.5-py36ha569be9_0 ...
    installing: chardet-3.0.4-py36h96c241c_1 ...
    installing: idna-2.6-py36h8628d0a_1 ...
    installing: pycosat-0.6.3-py36hee92d8f_0 ...
    installing: pycparser-2.18-py36h724b2fc_1 ...
    installing: pysocks-1.6.7-py36hfa33cec_1 ...
    installing: python.app-2-py36h54569d5_7 ...
    installing: ruamel_yaml-0.11.14-py36h9d7ade0_2 ...
    installing: six-1.11.0-py36h0e22d5e_1 ...
    installing: cffi-1.11.2-py36hd3e6348_0 ...
    installing: setuptools-36.5.0-py36h2134326_0 ...
    installing: cryptography-2.1.4-py36h842514c_0 ...
    installing: wheel-0.30.0-py36h5eb2c71_1 ...
    installing: pip-9.0.1-py36h1555ced_4 ...
    installing: pyopenssl-17.5.0-py36h51e4350_0 ...
    installing: urllib3-1.22-py36h68b9469_0 ...
    installing: requests-2.18.4-py36h4516966_1 ...
    installing: conda-4.3.31-py36_0 ...
    installation finished.
    Fetching package metadata .................
    Solving package specifications: .

    Package plan for installation in environment /Users/john/workspace/planemo/project_templates/seqtk_complete_cwl/cwltool_deps/_conda:

    The following packages will be UPDATED:

        conda: 4.3.31-py36_0 --> 4.3.33-py36_0 conda-forge

    conda-4.3.33-p 100% |#################################################################| Time: 0:00:00   1.13 MB/s


    Package plan for installation in environment /Users/john/workspace/planemo/project_templates/seqtk_complete_cwl/cwltool_deps/_conda/envs/__seqtk@1.2:

    The following NEW packages will be INSTALLED:

        seqtk: 1.2-1    bioconda
        zlib:  1.2.11-0 conda-forge


    [job seqtk_seq.cwl] completed success
    {
        "output1": {
            "checksum": "sha1$322e001e5a99f19abdce9f02ad0f02a17b5066c2",
            "basename": "out",
            "location": "file:///Users/john/workspace/planemo/project_templates/seqtk_complete_cwl/out",
            "path": "/Users/john/workspace/planemo/project_templates/seqtk_complete_cwl/out",
            "class": "File",
            "size": 150
        }
    }
    Final process status is success

This demonstrates that cwltool will install the packages needed on the first run, if we rerun cwltool it will
reuse that previous environment.

::

    $ cwltool --no-container --beta-conda-dependencies seqtk_seq.cwl seqtk_seq_job.yml
    /Users/john/workspace/planemo/.venv/bin/cwltool 1.0.20180508202931
    Resolved 'seqtk_seq.cwl' to 'file:///Users/john/workspace/planemo/project_templates/seqtk_complete_cwl/seqtk_seq.cwl'
    No handlers could be found for logger "rdflib.term"
    [job seqtk_seq.cwl] /private/tmp/docker_tmp4vvE_i$ seqtk \
        seq \
        -a \
        /private/var/folders/78/zxz5mz4d0jn53xf0l06j7ppc0000gp/T/tmpcvQ3Ph/stg2ef3a21c-9fb0-4099-88c2-36e24719901d/2.fastq > /private/tmp/docker_tmp4vvE_i/out
    [job seqtk_seq.cwl] completed success
    {
        "output1": {
            "checksum": "sha1$322e001e5a99f19abdce9f02ad0f02a17b5066c2",
            "basename": "out",
            "location": "file:///Users/john/workspace/planemo/project_templates/seqtk_complete_cwl/out",
            "path": "/Users/john/workspace/planemo/project_templates/seqtk_complete_cwl/out",
            "class": "File",
            "size": 150
        }
    }
    Final process status is success

And the same thing is possible with Toil_.

::

    $ cwltoil --no-container --beta-conda-dependencies seqtk_seq.cwl seqtk_seq_job.yml
    jlaptop17.local 2018-05-23 15:27:25,754 MainThread INFO toil.lib.bioio: Root logger is at level 'INFO', 'toil' logger at level 'INFO'.
    jlaptop17.local 2018-05-23 15:27:25,785 MainThread INFO toil.jobStores.abstractJobStore: The workflow ID is: '92328fb2-33b7-44cd-879f-41d8cbf94555'
    Resolved 'seqtk_seq.cwl' to 'file:///Users/john/workspace/planemo/project_templates/seqtk_complete_cwl/seqtk_seq.cwl'
    jlaptop17.local 2018-05-23 15:27:25,787 MainThread INFO cwltool: Resolved 'seqtk_seq.cwl' to 'file:///Users/john/workspace/planemo/project_templates/seqtk_complete_cwl/seqtk_seq.cwl'
    jlaptop17.local 2018-05-23 15:27:27,002 MainThread WARNING rdflib.term: http://schema.org/docs/!DOCTYPE html does not look like a valid URI, trying to serialize this will break.
    jlaptop17.local 2018-05-23 15:27:27,396 MainThread INFO rdflib.plugins.parsers.pyRdfa: Current options:
            preserve space                         : True
            output processor graph                 : True
            output default graph                   : True
            host language                          : RDFa Core
            accept embedded RDF                    : False
            check rdfa lite                        : False
            cache vocabulary graphs                : False

    jlaptop17.local 2018-05-23 15:27:29,797 MainThread INFO toil.common: Using the single machine batch system
    jlaptop17.local 2018-05-23 15:27:29,798 MainThread WARNING toil.batchSystems.singleMachine: Limiting maxCores to CPU count of system (8).
    jlaptop17.local 2018-05-23 15:27:29,798 MainThread WARNING toil.batchSystems.singleMachine: Limiting maxMemory to physically available memory (17179869184).
    jlaptop17.local 2018-05-23 15:27:29,808 MainThread INFO toil.common: Created the workflow directory at /var/folders/78/zxz5mz4d0jn53xf0l06j7ppc0000gp/T/toil-92328fb2-33b7-44cd-879f-41d8cbf94555-132281828025877
    jlaptop17.local 2018-05-23 15:27:29,808 MainThread WARNING toil.batchSystems.singleMachine: Limiting maxDisk to physically available disk (202669449216).
    jlaptop17.local 2018-05-23 15:27:29,815 MainThread INFO toil.common: User script ModuleDescriptor(dirPath='/Users/john/workspace/planemo/.venv/lib/python2.7/site-packages', name='toil.cwl.cwltoil', fromVirtualEnv=True) belongs to Toil. No need to auto-deploy it.
    jlaptop17.local 2018-05-23 15:27:29,816 MainThread INFO toil.common: No user script to auto-deploy.
    jlaptop17.local 2018-05-23 15:27:29,816 MainThread INFO toil.common: Written the environment for the jobs to the environment file
    jlaptop17.local 2018-05-23 15:27:29,816 MainThread INFO toil.common: Caching all jobs in job store
    jlaptop17.local 2018-05-23 15:27:29,816 MainThread INFO toil.common: 0 jobs downloaded.
    jlaptop17.local 2018-05-23 15:27:29,911 MainThread INFO toil: Running Toil version 3.15.0-0e3a87e738f5e0e7cff64bfdad337d592bd92704.
    jlaptop17.local 2018-05-23 15:27:29,911 MainThread INFO toil.realtimeLogger: Real-time logging disabled
    jlaptop17.local 2018-05-23 15:27:29,937 MainThread INFO toil.toilState: (Re)building internal scheduler state
    2018-05-23 15:27:29,937 - toil.toilState - INFO - (Re)building internal scheduler state
    jlaptop17.local 2018-05-23 15:27:29,938 MainThread INFO toil.leader: Found 1 jobs to start and 0 jobs with successors to run
    2018-05-23 15:27:29,938 - toil.leader - INFO - Found 1 jobs to start and 0 jobs with successors to run
    jlaptop17.local 2018-05-23 15:27:29,938 MainThread INFO toil.leader: Checked batch system has no running jobs and no updated jobs
    2018-05-23 15:27:29,938 - toil.leader - INFO - Checked batch system has no running jobs and no updated jobs
    jlaptop17.local 2018-05-23 15:27:29,938 MainThread INFO toil.leader: Starting the main loop
    2018-05-23 15:27:29,938 - toil.leader - INFO - Starting the main loop
    jlaptop17.local 2018-05-23 15:27:29,939 MainThread INFO toil.leader: Issued job 'file:///Users/john/workspace/planemo/project_templates/seqtk_complete_cwl/seqtk_seq.cwl' seqtk seq e/V/jobsxUpTU with job batch system ID: 0 and cores: 1, disk: 3.0 G, and memory: 2.0 G
    2018-05-23 15:27:29,939 - toil.leader - INFO - Issued job 'file:///Users/john/workspace/planemo/project_templates/seqtk_complete_cwl/seqtk_seq.cwl' seqtk seq e/V/jobsxUpTU with job batch system ID: 0 and cores: 1, disk: 3.0 G, and memory: 2.0 G
    jlaptop17.local 2018-05-23 15:27:31,409 MainThread INFO toil.leader: Job ended successfully: 'file:///Users/john/workspace/planemo/project_templates/seqtk_complete_cwl/seqtk_seq.cwl' seqtk seq e/V/jobsxUpTU
    2018-05-23 15:27:31,409 - toil.leader - INFO - Job ended successfully: 'file:///Users/john/workspace/planemo/project_templates/seqtk_complete_cwl/seqtk_seq.cwl' seqtk seq e/V/jobsxUpTU
    jlaptop17.local 2018-05-23 15:27:31,411 MainThread INFO toil.leader: Finished the main loop: no jobs left to run
    2018-05-23 15:27:31,411 - toil.leader - INFO - Finished the main loop: no jobs left to run
    jlaptop17.local 2018-05-23 15:27:31,411 MainThread INFO toil.serviceManager: Waiting for service manager thread to finish ...
    2018-05-23 15:27:31,411 - toil.serviceManager - INFO - Waiting for service manager thread to finish ...
    jlaptop17.local 2018-05-23 15:27:31,946 MainThread INFO toil.serviceManager: ... finished shutting down the service manager. Took 0.535056114197 seconds
    2018-05-23 15:27:31,946 - toil.serviceManager - INFO - ... finished shutting down the service manager. Took 0.535056114197 seconds
    jlaptop17.local 2018-05-23 15:27:31,947 MainThread INFO toil.statsAndLogging: Waiting for stats and logging collator thread to finish ...
    2018-05-23 15:27:31,947 - toil.statsAndLogging - INFO - Waiting for stats and logging collator thread to finish ...
    jlaptop17.local 2018-05-23 15:27:31,960 MainThread INFO toil.statsAndLogging: ... finished collating stats and logs. Took 0.0131621360779 seconds
    2018-05-23 15:27:31,960 - toil.statsAndLogging - INFO - ... finished collating stats and logs. Took 0.0131621360779 seconds
    jlaptop17.local 2018-05-23 15:27:31,961 MainThread INFO toil.leader: Finished toil run successfully
    2018-05-23 15:27:31,961 - toil.leader - INFO - Finished toil run successfully
    {
        "output1": {
            "checksum": "sha1$322e001e5a99f19abdce9f02ad0f02a17b5066c2",
            "basename": "out",
            "nameext": "",
            "nameroot": "out",
            "http://commonwl.org/cwltool#generation": 0,
            "location": "file:///Users/john/workspace/planemo/project_templates/seqtk_complete_cwl/out",
            "class": "File",
            "size": 150
        }
    jlaptop17.local 2018-05-23 15:27:31,972 MainThread INFO toil.common: Successfully deleted the job store: <toil.jobStores.fileJobStore.FileJobStore object at 0x10554d490>
    }2018-05-23 15:27:31,972 - toil.common - INFO - Successfully deleted the job store: <toil.jobStores.fileJobStore.FileJobStore object at 0x10554d490>

.. include:: _writing_conda_search.rst

----------------------------------------------------------------
Exercise - Leveraging Bioconda
----------------------------------------------------------------

Use the ``project_init`` command to download this exercise.

::

    $ planemo project_init --template conda_exercises_cwl conda_exercises
    $ cd conda_exercises/exercise1
    $ ls 
    pear.cwl              test-data

This project template contains a few exercises. The first uses a CWL tool for
`PEAR - Paired-End reAd mergeR <http://sco.h-its.org/exelixis/web/software/pear/>`__.
This tool however has no ``SoftwareRequirement`` or container annotations and so will
not work properly without modification.

1. Run ``planemo test pear.cwl`` to verify the tool does not function
   without dependencies defined.
2. Use ``--conda_requirements`` flag with ``planemo lint`` to verify it does
   indeed lack requirements.
3. Use ``planemo conda_search`` or the Anaconda_ website to search for the
   correct package and version in a best practice channel.
4. Update ``pear.cwl`` with the correct ``SoftwareRequirement`` hints.
5. Re-run the ``lint`` command from above to verify the tool now has the
   correct dependency definition.
6. Re-run the ``test`` command from above to verify the tool test now
   works properly.


.. include:: _writing_conda_new.rst

::

    $ planemo project_init --template conda_exercises_cwl conda_exercises
    $ cd conda_exercises/exercise2
    $ ls 
    fleeqtk_seq.cwl      fleeqtk_seq_tests.yml         test-data

.. include:: _writing_conda_fleeqtk.rst

1. Clone and branch Bioconda_.
2. Build a recipe for fleeqtk version 1.3. You may wish to use ``conda skeleton``, start from
   scratch, or copy the recipe of seqtk and work from there - any of these strategies should work.  
3. Use ``conda build`` or Bioconda tooling to build the recipe.
4. Run ``planemo conda_install --conda_use_local fleeqtk_seq.cwl`` to verify the resulting package
   can be built into a Galaxy environment.
5. Run ``planemo test fleeqtk_seq.cwl`` to verify the resulting package works as expected.

.. note: The planemo flag ``--conda_use_local`` causes Planemo to use locally built
     packages during dependency resolution and related commands.

.. include:: _writing_conda_recipe_complete.rst

.. _Software Requirements: https://www.commonwl.org/v1.0/CommandLineTool.html#SoftwareRequirement
.. _BioContainers: http://biocontainers.pro/
.. _Docker: https://www.docker.com/
.. _Singularity: https://singularity.lbl.gov/
.. _Common Workflow Language: https://www.commonwl.org/
.. _seqtk: https://github.com/lh3/seqtk
.. _fleeqtk: https://github.com/jmchilton/fleeqtk
.. _Bioconda: https://github.com/bioconda/bioconda-recipes
.. _Conda: https://conda.io/docs/
.. _Anaconda: https://anaconda.org/
.. _cwltool: https://github.com/common-workflow-language/cwltool
.. _Toil: https://github.com/BD2KGenomics/toil
.. _introductory tutorial: http://planemo.readthedocs.io/en/latest/writing_cwl.html
