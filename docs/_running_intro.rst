The Basics
================================

This tutorial will demonstrate workflow execution with a very simple test
workflow from the `workflow-testing <https://github.com/usegalaxy-eu/workflow-testing>`__
repository. This repository contains a number of workflows which are tested
regularly against the European Galaxy server.

::

    $ git clone https://github.com/usegalaxy-eu/workflow-testing.git
    $ cd workflow-testing/example3
    $ ls
    data  tutorial.ga  tutorial-job.yml  tutorial-tests.yml

The ``example3`` directory contains three files. Firstly, ``tutorial.ga``
contains a complete definition of the workflow in JSON format, which can be
easily exported from any Galaxy server after creating a new workflow. Secondly,
``tutorial-job.yml`` contains a list of the files and parameters (in YAML
format) which should be used for each of the workflow inputs upon execution.
Thirdly, ``tutorial-tests.yml`` contains a list of tests (similar to Galaxy
tool tests) which can be used to validate the outputs of the workflow.

The ``tutorial.ga`` workflow takes two input datasets and one input parameter,
and consists of two steps; firstly, ``Dataset 1`` and ``Dataset 2`` are
concatenated together, and secondly, a certain number of lines (specified by
the ``Number of lines`` parameter) are randomly selected. If you want to view
it in the Galaxy interface, you can do so with the command
``planemo workflow_edit tutorial.ga``.

The simplest way to run a workflow with planemo is on a locally hosted Galaxy
instance, just like executing a tool test with ``planemo test``. This can be
achieved with the command

::

    $ planemo run tutorial.ga tutorial-job.yml --output_directory . --output_json output.json


You can optionally (and probably should) add the ``--galaxy_root`` flag with
the location of a local copy of the Galaxy source code, which will allow the
instance to be spun up considerably faster.

Note that the ``--output_directory`` and ``--output_json`` flags are optional,
but allow saving the output to a local file. The contents should be something
like:

::

    $ cat tutorial_output.txt
    is
    hello
    world
    $ cat output.json
    {"output": {"class": "File", "path": "/home/user/workflow-testing/example3/tutorial_output.txt", "checksum": "sha1$4d7ab2b2bb0102ee5ec472a5971ca86081ff700c", "size": 15, "basename": "tutorial_output.txt", "nameroot": "tutorial_output", "nameext": ".txt"}}


You can also run the workflow on a local Dockerized Galaxy. For this, exactly
the same command can be used, with ``--engine docker_galaxy --ignore_dependency_problems``
appended. Please note that you need to have Docker installed and that it may take
significantly longer to complete than the previous command.

::

    $ planemo run tutorial.ga tutorial-job.yml --output_directory . --output_json output.json --engine docker_galaxy --ignore_dependency_problems


This introduces the concept of an engine, which Planemo provides to allow
workflows to be flexibly executed using the setup and workflow execution system
of the user's choice. The full list of engines provided by Galaxy is:
``galaxy`` (the default, used in the first example above), ``docker_galaxy``,
``cwltool``, ``toil`` and ``external_galaxy``.

As a final example to demonstrate workflow testing, try:

::

    $ planemo test tutorial.ga


In this case, Planemo automatically detects that it should test the workflow with
the ``tutorial-tests.yml``, so this file should be present and named correctly.
If you inspect its contents:

::

    $ cat tutorial-tests.yml
    - doc: Test outline for tutorial.ga
      job:
        Dataset 1:
          class: File
          path: "data/dataset1.txt"
        Dataset 2:
          class: File
          path: "data/dataset2.txt"
        Number of lines: 3
      outputs:
        output:
          class: File
          path: "data/output.txt"


you see that the job parameters are defined identically to the ``tutorial-job.yml``
file, with the addition of an output. For the test to pass, the output file
produced by the workflow must be identical to that stored in ``data/output.txt``.

The three commands above demonstrate the basics of workflow execution with
Planemo. For large scale workflow execution, however, it's likely that you would
prefer to use the more extensive resources provided by a public Galaxy server,
rather than running on a local instance. The tutorial therefore now turns to the
use of the ``galaxy_external`` engine, which as the name suggests, runs
workflows on a Galaxy external to Planemo.
