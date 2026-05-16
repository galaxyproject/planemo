====================================
Interacting with Galaxy workflows
====================================

Galaxy workflows are made up of a number of individual tools, which are
executed in sequence, automatically. They allow Galaxy users to perform complex
analyses made up of multiple simple steps.

Workflows can be easily created, edited and run using the Galaxy user interface
(i.e. in the web-browser), as is described in the
`workflow tutorial <https://training.galaxyproject.org/training-material/topics/galaxy-interface/tutorials/workflow-editor/tutorial.html>`__
provided by the Galaxy Training Network.

Planemo commands for interacting with workflows
===============================================

Planemo offers a number of convenient commands for interacting with Galaxy
workflows from the command line.
Here, you will use the following ones to interact with a small example workflow:

- ``planemo run``, which can be used to execute a Galaxy workflow with input datasets / parameters defined in a so-called *job file* in YAML format.

  If you are looking for a way to run workflows a very large number of times,
  or to automatically trigger workflow execution at particular times or as new
  data becomes available, this command is a great starting point!

- ``planemo test`` which cannot only be used to `test Galaxy tools <https://planemo.readthedocs.io/en/latest/writing_advanced.html#test-driven-development>`__, but also workflows.

  Similar to ``planemo run``, this command can be used to execute a Galaxy workflow,
  but it will also evaluate the success of the workflow execution by comparing workflow output datasets to expected results.

  This command enables test-driven development of Galaxy workflows.
  It can also form the basis of automated monitoring systems that, for example,
  check for compatibility between workflow versions and Galaxy server versions and instances.

  Input datasets / parameters and output assumptions are passed to this command in a *test file* in YAML format,
  which extends the job file format used with ``planemo run``.

- ``planemo workflow_job_init`` and ``planemo workflow_test_init``

  These are useful helper commands that generate templates of the *job file*
  expected by ``planemo run`` and of the *test file* expected by ``planemo test``,
  respectively, from a workflow definition file.

- ``planemo workflow_lint``, which lets you check a workflow for syntax errors and violation of workflow `best practices <https://planemo.readthedocs.io/en/latest/best_practices_workflows.html>`__.

- ``planemo list_invocations`` and ``planemo rerun``, which are great companions of ``planemo run``.

  ``planemo list_invocations`` provides information about the status of previous runs of a given workflow,
  while ``planemo rerun``, through its ``--invocation`` option lets you rerun failed jobs
  that resulted from any particular previous run of your workflow.

.. include:: _running_intro.rst
.. include:: _running_external.rst
