====================================
Running Galaxy workflows
====================================

Planemo offers a number of convenient commands for working with Galaxy
workflows. Workflows are made up of a number of individual tools, which are
executed in sequence, automatically. They allow Galaxy users to perform complex
analyses made up of multiple simple steps.

Workflows can be easily created, edited and run using the Galaxy user interface
(i.e. in the web-browser), as is described in the
`workflow tutorial <https://training.galaxyproject.org/training-material/topics/galaxy-interface/tutorials/workflow-editor/tutorial.html>`__
provided by the Galaxy Training Network. However, in some circumstances,
executing workflows may be awkward via the graphical interface. For example,
you might want to run workflows a very large number of times, or you might
want to automatically trigger workflow execution as a particular time as new
data becomes available. For these applications, being able to execute workflows
via the command line is very useful. This tutorial provides an introduction to
the ``planemo run`` command, which allows Galaxy tools and workflows to be
executed simply via the command line.

.. include:: _running_intro.rst
.. include:: _running_external.rst