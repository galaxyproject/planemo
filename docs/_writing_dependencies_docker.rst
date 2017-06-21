Dependencies and Docker
===========================================

For years Galaxy has supported running tools inside containers. The details
of how to run Galaxy tools inside of containers varies depending on the
Galaxy job runner but details can be found in Galaxy's job_conf.xml sample file.

This document doesn't describe how to run the containers, it describes how Galaxy
figures out which container to run for a given tool. There are currently
two strategies for finding containers for a tool - and they are each
discussed in detail in this document. The newer approach is more experimental
but will ultimately be considered the best practice approach - it is
to allow Galaxy to find or build a BioContainers_ container using ``requirement``
tags that resolve to best-practice Conda channels. The older approach is
to explicitly declare a container identifier in the tool XML.

While not as flexible as resolving arbitrary image IDs from URLs, the newer
approach has a few key advantages that make them a best practice:

- They provide superior reproducibility across Galaxy instances because the same
  binary Conda packages will automatically be used for both bare metal dependencies
  and inside containers.
- They are constructed automatically from existing Conda packages so tool
  developers shouldn't need to write ``Dockerfile`` s or register projects
  on Docker Hub.
- They are produced using mulled_ which produce very small containers
  that make deployment easy.

----------------------------------------------------------------
BioContainers_
----------------------------------------------------------------

.. note:: This section is a continuation of :ref:`dependencies_and_conda`,
    please review that section for background information on resolving
    requirements with Conda.

Finding and Building BioContainers_
----------------------------------------------------------------


If a tool contains requirements in best practice Conda channels, a
BioContainers_-style container can be found or built for it.

As reminder, ``planemo lint --conda_requirements <tool.xml>`` can be used
to check if a tool contains only best-practice ``requirement`` tags. The ``lint``
command can also be fed the ``--biocontainer`` flag to check if a
BioContainers_ container has been registered that is compatible with that tool.

Below is an example of using this with the completed ``seqtk_seq.xml``
tool from the introductory tutorial.

::

    $ planemo lint --biocontainers seqtk_seq.xml
    Linting tool /home/planemo/workspace/planemo/project_templates/seqtk_complete/seqtk_seq.xml
    Applying linter tests... CHECK
    .. CHECK: 1 test(s) found.
    Applying linter output... CHECK
    .. INFO: 1 outputs found.
    Applying linter inputs... CHECK
    .. INFO: Found 9 input parameters.
    Applying linter help... CHECK
    .. CHECK: Tool contains help section.
    .. CHECK: Help contains valid reStructuredText.
    Applying linter general... CHECK
    .. CHECK: Tool defines a version [0.1.0].
    .. CHECK: Tool defines a name [Convert to FASTA (seqtk)].
    .. CHECK: Tool defines an id [seqtk_seq].
    .. CHECK: Tool targets 16.01 Galaxy profile.
    Applying linter command... CHECK
    .. INFO: Tool contains a command.
    Applying linter citations... CHECK
    .. CHECK: Found 1 likely valid citations.
    Applying linter tool_xsd... CHECK
    .. INFO: File validates against XML schema.
    Applying linter biocontainer_registered... CHECK
    .. INFO: BioContainer best-practice container found [quay.io/biocontainers/seqtk:1.2--0].

This last line indicates that indeed a container has been registered
that is compatible with this tool -- ``quay.io/biocontainers/seqtk:1.2--0``.
We didn't do any extra work to build this container for this tool, all
BioConda_ recipes are packaged into containers and registered on quay.io_
as part of the BioContainers_ project.

This tool can be tested using ``planemo test`` in its BioContainer
Docker container using the flag ``--biocontainers`` as shown (in part) next.

::

    $ planemo test --biocontainers seqtk_seq.xml
    ...
    2017-03-01 08:18:19,669 INFO  [galaxy.tools.actions] Handled output named output1 for tool seqtk_seq (22.145 ms)
    2017-03-01 08:18:19,683 INFO  [galaxy.tools.actions] Added output datasets to history (14.604 ms)
    2017-03-01 08:18:19,703 INFO  [galaxy.tools.actions] Verified access to datasets for Job[unflushed,tool_id=seqtk_seq] (8.687 ms)
    2017-03-01 08:18:19,704 INFO  [galaxy.tools.actions] Setup for job Job[unflushed,tool_id=seqtk_seq] complete, ready to flush (20.380 ms)
    2017-03-01 08:18:19,719 INFO  [galaxy.tools.actions] Flushed transaction for job Job[id=2,tool_id=seqtk_seq] (15.191 ms)
    2017-03-01 08:18:20,120 INFO  [galaxy.jobs.handler] (2) Job dispatched
    2017-03-01 08:18:20,311 DEBUG [galaxy.tools.deps] Using dependency seqtk version 1.2 of type conda
    2017-03-01 08:18:20,312 DEBUG [galaxy.tools.deps] Using dependency seqtk version 1.2 of type conda
    2017-03-01 08:18:20,325 INFO  [galaxy.tools.deps.containers] Checking with container resolver [ExplicitDockerContainerResolver[]] found description [None]
    2017-03-01 08:18:20,468 INFO  [galaxy.tools.deps.containers] Checking with container resolver [CachedMulledDockerContainerResolver[namespace=None]] found description [None]
    2017-03-01 08:18:20,881 INFO  [galaxy.tools.deps.containers] Checking with container resolver [MulledDockerContainerResolver[namespace=biocontainers]] found description [ContainerDescription[identifier=quay.io/biocontainers/seqtk:1.2--0,type=docker]]
    2017-03-01 08:18:20,904 INFO  [galaxy.jobs.command_factory] Built script [/tmp/tmpw8_UQm/job_working_directory/000/2/tool_script.sh] for tool command [seqtk seq -a '/tmp/tmpw8_UQm/files/000/dataset_1.dat' > '/tmp/tmpw8_UQm/files/000/dataset_2.dat']
    2017-03-01 08:18:21,060 DEBUG [galaxy.tools.deps] Using dependency samtools version None of type conda
    2017-03-01 08:18:21,061 DEBUG [galaxy.tools.deps] Using dependency samtools version None of type conda
    ok
    
    ----------------------------------------------------------------------
    XML: /private/tmp/tmpw8_UQm/xunit.xml
    ----------------------------------------------------------------------
    Ran 1 test in 11.926s
    
    OK
    2017-03-01 08:18:26,726 INFO  [test_driver] Shutting down
    ...
    2017-03-01 08:18:26,732 INFO  [galaxy.jobs.handler] job handler stop queue stopped
    Testing complete. HTML report is in "/home/planemo/workspace/planemo/tool_test_output.html".
    All 1 test(s) executed passed.
    seqtk_seq[0]: passed
    $

A very important line here is::

    2017-03-01 08:18:20,881 INFO  [galaxy.tools.deps.containers] Checking with container resolver [MulledDockerContainerResolver[namespace=biocontainers]] found description [ContainerDescription[identifier=quay.io/biocontainers/seqtk:1.2--0,type=docker]]

This line indicates that Galaxy was able to find a container for this tool and
executed the tool in that container.

For interactive testing, the ``planemo serve`` command also implements the
``--biocontainers`` flag.

In this seqtk example the relevant BioContainer already existed on quay.io_,
this won't always be the case. For tools that contain multiple ``requirement``
tags an existing container likely won't exist. The mulled_ toolkit
(distributed with planemo or available standalone) can be used to build
containers for such tools. For such tools, if Galaxy is configured to use
BioContainers it will attempt to build these containers on the fly by default
(though this behavior can be disabled).

You can try it directly using the ``mull`` command in Planemo. The ``conda_testing``
Planemo project template has a toy example tool with two requirements for
demonstrating this - `bwa_and_samtools.xml
<https://github.com/galaxyproject/planemo/blob/master/project_templates/conda_testing/bwa_and_samtools.xml>`__.

::

    $ planemo project_init --template=conda_testing conda_testing
    $ cd conda_testing/
    $ planemo mull bwa_and_samtools.xml
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

As the output indicates, this command built the container named
``quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:03dc1d2818d9de56938078b8b78b82d967c1f820-0``.
This is the same namespace / URL that would be used if or when published by
the BioContainers_ project.

.. note:: The first part of this ``mulled-v2`` hash is a hash of the package names
    that went into it, the second the packages used and build number. Check out
    the `Multi-package Containers <http://biocontainers.pro/multi-package-containers/>`__
    web application to explore best practice channels and build such hashes.

We can see this new container when running the Docker command ``images`` and
explore the new container interactively with ``docker run``.

::

    $ docker images
    REPOSITORY                                                                 TAG                                          IMAGE ID            CREATED              SIZE
    quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40   03dc1d2818d9de56938078b8b78b82d967c1f820-0   a740fe1e6a9e        16 hours ago         104 MB
    quay.io/biocontainers/seqtk                                                1.2--0                                       10bc359ebd30        2 days ago           7.34 MB
    continuumio/miniconda                                                      latest                                       6965a4889098        3 weeks ago          437 MB
    bgruening/busybox-bash                                                     0.1                                          3d974f51245c        9 months ago         6.73 MB
    $ docker run -i -t quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:03dc1d2818d9de56938078b8b78b82d967c1f820-0 /bin/bash
    bash-4.2# which samtools
    /usr/local/bin/samtools
    bash-4.2# which bwa
    /usr/local/bin/bwa

As before, we can test running the tool inside its container in Galaxy using
the ``--biocontainers`` flag.

::

    $ planemo test --biocontainers bwa_and_samtools.xml
    ...
    2017-03-01 10:20:58,077 INFO  [galaxy.tools.actions] Handled output named output_2 for tool bwa_and_samtools (17.443 ms)
    2017-03-01 10:20:58,090 INFO  [galaxy.tools.actions] Added output datasets to history (12.935 ms)
    2017-03-01 10:20:58,095 INFO  [galaxy.tools.actions] Verified access to datasets for Job[unflushed,tool_id=bwa_and_samtools] (0.021 ms)
    2017-03-01 10:20:58,096 INFO  [galaxy.tools.actions] Setup for job Job[unflushed,tool_id=bwa_and_samtools] complete, ready to flush (5.755 ms)
    2017-03-01 10:20:58,116 INFO  [galaxy.tools.actions] Flushed transaction for job Job[id=1,tool_id=bwa_and_samtools] (19.582 ms)
    2017-03-01 10:20:58,869 INFO  [galaxy.jobs.handler] (1) Job dispatched
    2017-03-01 10:20:59,067 DEBUG [galaxy.tools.deps] Using dependency bwa version 0.7.15 of type conda
    2017-03-01 10:20:59,067 DEBUG [galaxy.tools.deps] Using dependency samtools version 1.3.1 of type conda
    2017-03-01 10:20:59,067 DEBUG [galaxy.tools.deps] Using dependency bwa version 0.7.15 of type conda
    2017-03-01 10:20:59,068 DEBUG [galaxy.tools.deps] Using dependency samtools version 1.3.1 of type conda
    2017-03-01 10:20:59,083 INFO  [galaxy.tools.deps.containers] Checking with container resolver [ExplicitContainerResolver[]] found description [None]
    2017-03-01 10:20:59,142 INFO  [galaxy.tools.deps.containers] Checking with container resolver [CachedMulledDockerContainerResolver[namespace=biocontainers]] found description [ContainerDescription[identifier=quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:03dc1d2818d9de56938078b8b78b82d967c1f820-0,type=docker]]
    2017-03-01 10:20:59,163 INFO  [galaxy.jobs.command_factory] Built script [/tmp/tmpQs0gyp/job_working_directory/000/1/tool_script.sh] for tool command [bwa > /tmp/tmpQs0gyp/files/000/dataset_1.dat 2>&1 ; samtools > /tmp/tmpQs0gyp/files/000/dataset_2.dat 2>&1]
    2017-03-01 10:20:59,367 DEBUG [galaxy.tools.deps] Using dependency samtools version None of type conda
    2017-03-01 10:20:59,367 DEBUG [galaxy.tools.deps] Using dependency samtools version None of type conda
    ok
    
    ----------------------------------------------------------------------
    XML: /private/tmp/tmpQs0gyp/xunit.xml
    ----------------------------------------------------------------------
    Ran 1 test in 7.553s
    
    OK
    2017-03-01 10:21:05,223 INFO  [test_driver] Shutting down
    2017-03-01 10:21:05,224 INFO  [test_driver] Shutting down embedded galaxy web server
    2017-03-01 10:21:05,226 INFO  [test_driver] Embedded web server galaxy stopped
    2017-03-01 10:21:05,226 INFO  [test_driver] Stopping application galaxy
    ...
    2017-03-01 10:21:05,228 INFO  [galaxy.jobs.handler] job handler stop queue stopped
    Testing complete. HTML report is in "/home/planemo/workspace/planemo/tool_test_output.html".
    All 1 test(s) executed passed.
    bwa_and_samtools[0]: passed

In particular take note of the line::

    2017-03-01 10:20:59,142 INFO  [galaxy.tools.deps.containers] Checking with container resolver [CachedMulledDockerContainerResolver[namespace=biocontainers]] found description [ContainerDescription[identifier=quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:03dc1d2818d9de56938078b8b78b82d967c1f820-0,type=docker]]

Here we can see the container ID (``quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:03dc1d2818d9de56938078b8b78b82d967c1f820-0``)
from earlier has been cached on our Docker host is picked up by Galaxy. This is used to run the simple
tool tests and indeed they pass.

In our initial seqtk example, the container resolver that matched was of type
``MulledDockerContainerResolver`` indicating that the Docker image would be downloaded
from the BioContainer repository and this time the resolve that matched was of type
``CachedMulledDockerContainerResolver`` meaning that Galaxy would just use the locally
cached version from the Docker host (i.e. the one we built with ``planemo mull``
above).

Planemo doesn't yet expose options that make it possible to build mulled
containers for local packages that have yet to be published to anaconda.org
but the mulled toolkit allows this. See mulled_ documentation for more
information. However, once a container for a local package is built with
``mulled-build-tool`` the ``--biocontainers`` command should work to test
it.


Publishing BioContainers_
----------------------------------------------------------------

Building unpublished BioContainers on the fly is great for testing but
for production use and to increase reproducibility such containers should
ideally be published as well.

BioContainers_ maintains a registry of package combinations to be published
using these long mulled hashes. This registry is represented as a Github repository
named `multi-package-containers <https://github.com/biocontainers/multi-package-containers>`__.
The Planemo command ``container_register`` will inspect a tool and open a
Github pull request to add the tool's combination
of packages to the registry. Once merged, this pull request will
result in the corresponding BioContainers image to be published (with the
correct mulled has as its name) - these can be subsequently be picked up by
Galaxy.

Various Github related settings need to be configured in order for Planemo
to be able to open pull requests on your behalf as part of the
``container_register`` command. To simplify all of this - the Planemo community
maintains a list of Github repositories containing Galaxy and/or CWL tools that
are scanned daily by Travis_. For each such repository, the Travis job will run
``container_register`` across the repository on all tools resulting in new registry
pull requests for all new combinations of tools. This list is maintained
in a script named ``monitor.sh`` in the `planemo-monitor
<https://github.com/galaxyproject/planemo-monitor/>`__ repository. The easiest way
to ensure new containers are built for your tools is simply to open open a pull
request to add your tool repositories to this list.

----------------------------------------------------------------
Explicit Annotation
----------------------------------------------------------------

This section of documentation needs to be filled out but a detailed
example is worked through `this documentation
<https://github.com/apetkau/galaxy-hackathon-2014>`__ from Aaron Petkau
(@apetkau) built at the 2014 Galaxy Community Conference Hackathon.

.. _BioContainers: http://biocontainers.pro/
.. _mulled: https://github.com/BioContainers/auto-mulled
.. _quay.io: https://quay.io
.. _Travis: https://travis-ci.org