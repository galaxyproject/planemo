The Basics
================================

This tutorial will demonstrate workflow execution with a very simple test
workflow from the `workflow-testing <https://github.com/usegalaxy-eu/workflow-testing>`__
repository. This repository contains a number of workflows which are tested
regularly against the European Galaxy server.

::

    $ git clone https://github.com/usegalaxy-eu/workflow-testing.git
    $ cd example1
    $ ls
    wf3-shed-tools.ga  wf3-shed-tools-job.yml  wf3-shed-tools-test.yml

The `example1` directory contains three files. Firstly, `wf-shed-tools.ga`
contains a complete definition of the workflow in JSON format, which can be
easily exported from any Galaxy server after creating a new workflow. Secondly,
`wf3-shed-tools-job.yml` contains a list of the files and parameters (in YAML
format) which should be used for each of the workflow inputs upon execution.
Thirdly, `wf3-shed-tools-test.yml` contains a list of tests (similar to Galaxy
tool tests) which can be used to validate the outputs of the workflow.

The simplest way to run a workflow with planemo is on a locally hosted Galaxy
instance, just like executing a tool test with `planemo test`. This can be
achieved with the command

::

    $ planemo run wf3-shed-tools.ga wf3-shed-tools-job.yml


You can optionally add the `--galaxy_root` flag with the location of a local
copy of the Galaxy source code, which should allow the instance to be spun up
considerably faster.

You can also run the workflow on a local Dockerized Galaxy. For this, exactly
the same command can be used, with `--engine docker_galaxy` appended:

::

    $ planemo run wf3-shed-tools.ga wf3-shed-tools-job.yml --engine docker_galaxy


This introduces the concept of an engine, which Planemo provides to allow
workflows to be flexibly executed using the workflow management system or setup
of the user's choice. The full list of engines provided by Galaxy is:
`galaxy` (the default, used in the first example above), `docker_galaxy`,
`cwltool`, `toil` and `external_galaxy`.

As a final example to demonstrate workflow testing, try:

::

    $ planemo test wf3-shed-tools.ga wf3-shed-tools-job.yml --engine docker_galaxy



The three commands above demonstrate the basics of workflow execution with
Planemo. For large scale workflow execution, however, it's likely that you would
prefer to use the more extensive resources provided by a public Galaxy server,
rather than running on a local instance. The tutorial therefore now turns to the
use of the `galaxy_external` engine, which as the name suggests, uses runs
workflows on a Galaxy external to Planemo.