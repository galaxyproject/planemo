.. _dependencies_and_containers_cwl:

Dependencies and Containers
===========================================

.. note:: This section is a continuation of :ref:`dependencies_and_conda_cwl`,
    please review that section for background information on resolving
    `Software Requirements`_ with Conda.

Common Workflow Language tools can be annotated with arbitrary Docker requirements,
see the `CWL User Guide <http://www.commonwl.org/user_guide/07-containers/>`__
for a discussion about how to do this in general.

This document will discuss some techniques to find containers automatically from
the ``SoftwareRequirement`` annotations when using Planemo, cwltool_, or Toil_.
You will ultimately want to explicitly annotate your tools with the containers
we describe here so that other CWL implementations will be able to find containers
for your tool, but there are real advantages to using these containers instead
of ad-hoc things you may build with a ``Dockerfile``.

- They provide superior reproducibility because the same binary Conda packages
  will automatically be used for both bare metal dependencies and inside containers.
- They are constructed automatically from existing Conda packages so you as a tool
  developer won't need to write ``Dockerfile`` s or register projects on Docker Hub.
- They are produced using mulled_ which produce very small containers
  that make deployment easier regardless of the CWL implementation you are using.
- Annotating `Software Requirements`_ reduces the opaqueness of the Docker process.
  With this method it is entirely traceable how the container was constructed from
  what sources were fetched, which exact build of every dependency was used, to how
  packages in the container were built. Beyond that metadata about the packages can be
  fetched from Bioconda_ (e.g. `this
  <https://github.com/BioContainers/biotools-bioconda-ids/blob/master/mapping.csv>`__).

Read more about this reproducibility stack in our preprint `Practical computational
reproducibility in the life sciences <https://www.biorxiv.org/content/early/2017/10/10/200683>`__.

----------------------------------------------------------------
BioContainers_
----------------------------------------------------------------

.. note:: This section is a continuation of :ref:`dependencies_and_conda_cwl`,
    please review that section for background information on resolving
    `Software Requirements`_ with Conda.

Finding BioContainers_
----------------------------------------------------------------

If a tool contains `Software Requirements`_ in best practice Conda channels, a
BioContainers_-style container can be found or built for it.

As reminder, ``planemo lint --conda_requirements <tool.cwl>`` can be used
to check if a tool contains only best-practice ``requirement`` tags. The ``lint``
command can also be fed the ``--biocontainers`` flag to check if a
BioContainers_ container has been registered that is compatible with that tool.

The Conda exercises project template has an example tool (``exercise3``) that we
can use to demonstrate ``--biocontainers``. If you are continuing from the Conda
tutorial, simply move to ``../exercise3`` otherwise using ``planemo project_init``
to grab the exercise as show below.

::

    $ planemo project_init --template conda_exercises_cwl conda_exercises 
    $ cd conda_exercises/exercise3
    $ planemo lint --biocontainers seqtk_seq.cwl
    Linting tool /home/planemo/conda_exercises_cwl/exercise_3/seqtk_seq.cwl
    Applying linter general... CHECK
    .. CHECK: Tool defines a version [0.0.1].
    .. CHECK: Tool defines a name [Convert to FASTA (seqtk)].
    .. CHECK: Tool defines an id [seqtk_seq].
    .. CHECK: Tool specifies profile version [16.04].
    Applying linter cwl_validation... CHECK
    .. INFO: CWL appears to be valid.
    Applying linter docker_image... WARNING
    .. WARNING: Tool does not specify a DockerPull source.
    Applying linter new_draft... CHECK
    .. INFO: Modern CWL version [v1.0]
    Applying linter biocontainer_registered... CHECK
    .. INFO: BioContainer best-practice container found [quay.io/biocontainers/seqtk:1.2--1].
    Failed linting

.. include:: _writing_containers_linter_explain.rst

::

    $ planemo test --biocontainers seqtk_seq.cwl
    Enable beta testing mode for testing.
    cwltool INFO: /Users/john/workspace/planemo/.venv/bin/planemo 1.0.20180508202931
    cwltool INFO: Resolved '/Users/john/workspace/planemo/project_templates/conda_exercises_cwl/exercise_3/seqtk_seq.cwl' to 'file:///Users/john/workspace/planemo/project_templates/conda_exercises_cwl/exercise_3/seqtk_seq.cwl'
    galaxy.tools.deps.containers INFO: Checking with container resolver [ExplicitContainerResolver[]] found description [None]
    galaxy.tools.deps.containers INFO: Checking with container resolver [CachedMulledDockerContainerResolver[namespace=biocontainers]] found description [None]
    galaxy.tools.deps.containers INFO: Checking with container resolver [MulledDockerContainerResolver[namespace=biocontainers]] found description [ContainerDescription[identifier=quay.io/biocontainers/seqtk:1.2--1,type=docker]]
    cwltool INFO: [job seqtk_seq.cwl] /private/tmp/docker_tmpMEipaU$ docker \
        run \
        -i \
        --volume=/private/tmp/docker_tmpMEipaU:/private/tmp/docker_tmpMEipaU:rw \
        --volume=/private/var/folders/78/zxz5mz4d0jn53xf0l06j7ppc0000gp/T/tmpxkm9dp:/tmp:rw \
        --volume=/Users/john/workspace/planemo/project_templates/conda_exercises_cwl/exercise_3/test-data/2.fastq:/private/var/folders/78/zxz5mz4d0jn53xf0l06j7ppc0000gp/T/tmpjAVM_1/stgddf6fc2a-dd13-4322-9b88-68571a1697dd/2.fastq:ro \
        --workdir=/private/tmp/docker_tmpMEipaU \
        --read-only=true \
        --log-driver=none \
        --user=502:20 \
        --rm \
        --env=TMPDIR=/tmp \
        --env=HOME=/private/tmp/docker_tmpMEipaU \
        quay.io/biocontainers/seqtk:1.2--1 \
        seqtk \
        seq \
        -a \
        /private/var/folders/78/zxz5mz4d0jn53xf0l06j7ppc0000gp/T/tmpjAVM_1/stgddf6fc2a-dd13-4322-9b88-68571a1697dd/2.fastq > /private/tmp/docker_tmpMEipaU/out
    cwltool INFO: [job seqtk_seq.cwl] completed success
    cwltool INFO: Final process status is success
    All 1 test(s) executed passed.
    seqtk_seq_0: passed

----------------------------------------------------------------
Exercise - Leveraging Bioconda
----------------------------------------------------------------

1. Try the above command without the ``--biocontainers`` argument. Verify the tool does not
   run in a container by default.
2. Add a `DockerRequirement <http://www.commonwl.org/user_guide/07-containers/>`__ based on the
   the lint output above to annotate this tool with a Biocontainers Docker container and
   rerun test to verify the tool works now.

Building BioContainers_
----------------------------------------------------------------

In this seqtk example above the relevant BioContainer already existed on quay.io_,
this won't always be the case. For tools that contain multiple `Software Requirements`_
tags an existing container likely won't exist. The mulled_ toolkit
(distributed with planemo or available standalone) can be used to build
containers for such tools. For such tools, if cwltool_ or Toil_ is configured to use
BioContainers it will attempt to build these containers on the fly by default
(though this behavior can be disabled).

You can try it directly using the ``mull`` command in Planemo. The ``conda_testing``
Planemo project template has a toy example tool with two requirements for
demonstrating this - `bwa_and_samtools.cwl
<https://github.com/galaxyproject/planemo/blob/master/project_templates/conda_testing_cwl/bwa_and_samtools.cwl>`__.

::

    $ planemo project_init --template=conda_testing_cwl conda_testing
    $ cd conda_testing/
    $ planemo mull bwa_and_samtools.cwl
    /Users/john/.planemo/involucro -v=3 -f /Users/john/workspace/planemo/.venv/lib/python2.7/site-packages/galaxy_lib-17.9.0-py2.7.egg/galaxy/tools/deps/mulled/invfile.lua -set CHANNELS='iuc,bioconda,r,defaults,conda-forge' -set TEST='true' -set TARGETS='samtools=1.3.1,bwa=0.7.15' -set REPO='quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:03dc1d2818d9de56938078b8b78b82d967c1f820' -set BINDS='build/dist:/usr/local/' -set PREINSTALL='conda install --quiet --yes conda=4.3' build
    /Users/john/.planemo/involucro -v=3 -f /Users/john/workspace/planemo/.venv/lib/python2.7/site-packages/galaxy_lib-17.9.0-py2.7.egg/galaxy/tools/deps/mulled/invfile.lua -set CHANNELS='iuc,bioconda,r,defaults,conda-forge' -set TEST='true' -set TARGETS='samtools=1.3.1,bwa=0.7.15' -set REPO='quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:03dc1d2818d9de56938078b8b78b82d967c1f820' -set BINDS='build/dist:/usr/local/' -set PREINSTALL='conda install --quiet --yes conda=4.3' build
    [Jun 19 11:28:35] DEBU Run file [/Users/john/workspace/planemo/.venv/lib/python2.7/site-packages/galaxy_lib-17.9.0-py2.7.egg/galaxy/tools/deps/mulled/invfile.lua]
    [Jun 19 11:28:35] STEP Run image [continuumio/miniconda:latest] with command [[rm -rf /data/dist]]
    [Jun 19 11:28:35] DEBU Creating container [step-730a02d79e]
    [Jun 19 11:28:35] DEBU Created container [5e4b5f83c455 step-730a02d79e], starting it
    [Jun 19 11:28:35] DEBU Container [5e4b5f83c455 step-730a02d79e] started, waiting for completion
    [Jun 19 11:28:36] DEBU Container [5e4b5f83c455 step-730a02d79e] completed with exit code [0] as expected
    [Jun 19 11:28:36] DEBU Container [5e4b5f83c455 step-730a02d79e] removed
    [Jun 19 11:28:36] STEP Run image [continuumio/miniconda:latest] with command [[/bin/sh -c conda install --quiet --yes conda=4.3 && conda install  -c iuc -c bioconda -c r -c defaults -c conda-forge  samtools=1.3.1 bwa=0.7.15 -p /usr/local --copy --yes --quiet]]
    [Jun 19 11:28:36] DEBU Creating container [step-e95bf001c8]
    [Jun 19 11:28:36] DEBU Created container [72b9ca0e56f8 step-e95bf001c8], starting it
    [Jun 19 11:28:37] DEBU Container [72b9ca0e56f8 step-e95bf001c8] started, waiting for completion
    [Jun 19 11:28:46] SOUT Fetching package metadata .........
    [Jun 19 11:28:47] SOUT Solving package specifications: .
    [Jun 19 11:28:50] SOUT
    [Jun 19 11:28:50] SOUT Package plan for installation in environment /opt/conda:
    [Jun 19 11:28:50] SOUT
    [Jun 19 11:28:50] SOUT The following packages will be UPDATED:
    [Jun 19 11:28:50] SOUT
    [Jun 19 11:28:50] SOUT conda: 4.3.11-py27_0 --> 4.3.22-py27_0
    [Jun 19 11:28:50] SOUT
    [Jun 19 11:29:04] SOUT Fetching package metadata .................
    [Jun 19 11:29:06] SOUT Solving package specifications: .
    [Jun 19 11:29:56] SOUT
    [Jun 19 11:29:56] SOUT Package plan for installation in environment /usr/local:
    [Jun 19 11:29:56] SOUT
    [Jun 19 11:29:56] SOUT The following NEW packages will be INSTALLED:
    [Jun 19 11:29:56] SOUT
    [Jun 19 11:29:56] SOUT bwa:        0.7.15-1      bioconda
    [Jun 19 11:29:56] SOUT curl:       7.52.1-0
    [Jun 19 11:29:56] SOUT libgcc:     5.2.0-0
    [Jun 19 11:29:56] SOUT openssl:    1.0.2l-0
    [Jun 19 11:29:56] SOUT pip:        9.0.1-py27_1
    [Jun 19 11:29:56] SOUT python:     2.7.13-0
    [Jun 19 11:29:56] SOUT readline:   6.2-2
    [Jun 19 11:29:56] SOUT samtools:   1.3.1-5       bioconda
    [Jun 19 11:29:56] SOUT setuptools: 27.2.0-py27_0
    [Jun 19 11:29:56] SOUT sqlite:     3.13.0-0
    [Jun 19 11:29:56] SOUT tk:         8.5.18-0
    [Jun 19 11:29:56] SOUT wheel:      0.29.0-py27_0
    [Jun 19 11:29:56] SOUT zlib:       1.2.8-3
    [Jun 19 11:29:56] SOUT
    [Jun 19 11:29:57] DEBU Container [72b9ca0e56f8 step-e95bf001c8] completed with exit code [0] as expected
    [Jun 19 11:29:57] DEBU Container [72b9ca0e56f8 step-e95bf001c8] removed
    [Jun 19 11:29:57] STEP Wrap [build/dist] as [quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:03dc1d2818d9de56938078b8b78b82d967c1f820-0]
    [Jun 19 11:29:57] DEBU Creating container [step-6f1c176372]
    [Jun 19 11:29:58] DEBU Packing succeeded


.. include:: _writing_containers_mulled_output.rst

As before, we can test running the tool inside its container in cwltool_ using
the ``--biocontainers`` flag.

::

    $ planemo test --biocontainers bwa_and_samtools.cwl
    Enable beta testing mode for testing.
    cwltool INFO: /Users/john/workspace/planemo/.venv/bin/planemo 1.0.20180508202931
    cwltool INFO: Resolved '/Users/john/workspace/planemo/project_templates/conda_testing_cwl/bwa_and_samtools.cwl' to 'file:///Users/john/workspace/planemo/project_templates/conda_testing_cwl/bwa_and_samtools.cwl'
    galaxy.tools.deps.containers INFO: Checking with container resolver [ExplicitContainerResolver[]] found description [None]
    galaxy.tools.deps.containers INFO: Checking with container resolver [CachedMulledDockerContainerResolver[namespace=biocontainers]] found description [ContainerDescription[identifier=quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:03dc1d2818d9de56938078b8b78b82d967c1f820-0,type=docker]]
    cwltool INFO: [job bwa_and_samtools.cwl] /private/tmp/docker_tmpYJnmO4$ docker \
        run \
        -i \
        --volume=/private/tmp/docker_tmpYJnmO4:/private/tmp/docker_tmpYJnmO4:rw \
        --volume=/private/var/folders/78/zxz5mz4d0jn53xf0l06j7ppc0000gp/T/tmpVI06me:/tmp:rw \
        --workdir=/private/tmp/docker_tmpYJnmO4 \
        --read-only=true \
        --user=502:20 \
        --rm \
        --env=TMPDIR=/tmp \
        --env=HOME=/private/tmp/docker_tmpYJnmO4 \
        quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:03dc1d2818d9de56938078b8b78b82d967c1f820-0 \
        sh \
        -c \
        'bwa > bwa_help.txt 2>&1; samtools > samtools_help.txt 2>&1'
    cwltool INFO: [job bwa_and_samtools.cwl] completed success
    cwltool INFO: Final process status is success
    All 1 test(s) executed passed.
    bwa_and_samtools_0: passed

In particular take note of the line::

    2017-03-01 10:20:59,142 INFO  [galaxy.tools.deps.containers] Checking with container resolver [CachedMulledDockerContainerResolver[namespace=biocontainers]] found description [ContainerDescription[identifier=quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:03dc1d2818d9de56938078b8b78b82d967c1f820-0,type=docker]]

Here we can see the container ID (``quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:03dc1d2818d9de56938078b8b78b82d967c1f820-0``)
from earlier has been cached on our Docker host is picked up by cwltool_. This is used to run the
simple tool tests and indeed they pass.

In our initial seqtk_ example, the container resolver that matched was of type
``MulledDockerContainerResolver`` indicating that the Docker image would be downloaded
from the BioContainers_ repository and this time the resolve that matched was of type
``CachedMulledDockerContainerResolver`` meaning that cwltool_ would just use the locally
cached version from the Docker host (i.e. the one we built with ``planemo mull``
above).

.. include:: _writing_containers_mulled_build.rst

.. include:: _writing_containers_publishing.rst

.. _Software Requirements: https://www.commonwl.org/v1.0/CommandLineTool.html#SoftwareRequirement
.. _BioContainers: http://biocontainers.pro/
.. _mulled: https://github.com/BioContainers/auto-mulled
.. _quay.io: https://quay.io
.. _cwltool: https://github.com/common-workflow-language/cwltool
.. _Toil: https://github.com/BD2KGenomics/toil
