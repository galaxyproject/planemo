Automating Galaxy Workflows
=============================

The BioBlend_ library can be used to invoke and monitor Galaxy workflows.
Planemo provides a higher-level interfaces to working with workflows -
both as a Python library and via the command line. Planemo can take care of
orchestrating details such as launching and configuring a Galaxy instance,
installing required tools for the workflow, monitoring the workflow invocation,
and downloading the results back to a local directory.

::

    planemo run workflow.ga job.yml

The format of the job file should be YAML or JSON and matches the job
definition used by the workflow `test format
<https://planemo.readthedocs.io/en/latest/test_format.html>`__.

A job template for a particular workflow can be created with the
``workflow_job_init`` command.

::

    planemo workflow_job_init my-workflow.ga

This will create a ``my-workflow-job.yml`` file in the current directory.

Planemo run and the underlying Python code when used as a library support
a wide range of arguments for how to run a Galaxy workflow.

.. _BioBlend: https://github.com/galaxyproject/bioblend
