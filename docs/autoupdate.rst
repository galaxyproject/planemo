=============================
Autoupdating tools
=============================

Galaxy tools use underlying command-line software which are specified using ``<requirement>`` tags in the XML wrapper. If developers continue to release new versions of the software, the creator/maintainer needs to ensure that the tool requirements are bumped, otherwise the tool will become updated.

Planemo provides a ``autoupdate`` subcommand which can be used to perform this task automatically. Basic usage is as follows:

::

    planemo autoupdate tool_wrapper.xml

There are various flags which can be applied; these are some of the most important:
  - ``--recursive``, which performs the autoupdate recursively on subdirectories
  - ``--dry-run``, which checks if tool requirements are up-to-date without making the necessary change automatically
  - ``--test``, which runs tests on all autoupdated tools. If this flag is used, all options available for ``planemo test`` are also available.
  - ``--update_test_data`` (if ``--test`` is also selected) which will update test data for failed tests
  - ``--skiplist``, pointing to a list of tool wrappers for which updates should be skipped
  - ``--skip_requirements`` with a comma-separated list of packages not to update.  ``python``, ``perl``, ``r-base`` are skipped by default.

One of the most efficient ways to use it is by implementing a CI cron job which runs the command on an entire GitHub repository of tool wrappers.

Formatting tools
=============================

``autoupdate`` assumes that XML tools are formatted in a certain way - in accordance with the `IUC best practices`_. In particular:

  - the tool ``version`` attribute must be set to ``@TOOL_VERSION@+galaxy0`` or ``@TOOL_VERSION@+galaxy@VERSION_SUFFIX@``
  - A token ``@TOOL_VERSION@`` should be created which corresponds to the version number of the main requirement.
  - Optionally, a token ``@VERSION_SUFFIX@`` should be created, which should be an integer representing the number of times the XML wrapper has been updated since ``@TOOL_VERSION@`` was updated.

Updating workflows
=============================

The ``autoupdate`` subcommand can also be used to automatically update workflows so that they are using the most recent Galaxy tools available.

::

    planemo autoupdate workflow.ga

In the basic usage, a local Galaxy instance will be spun up and the workflow uploaded, refactored to include the most recent tool versions, and re-downloaded.

Workflows can also be updated against an external galaxy, for example:

::

    planemo autoupdate workflow.ga --profile usegalaxy-eu

In this case, the workflow returned will contain the most recent tool and subworkflow versions available on that Galaxy server.

Implementing an autoupdate CI job
=================================

 A useful way to use the autoupdate command is to implement it as a CI job, so that tools in a repo can be updated on a regular basis (e.g. weekly). An example implementation is the `planemo-autoupdate`_ GitHub bot.


.. _IUC best practices: https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/tool_xml.html
.. _planemo-autoupdate: https://github.com/planemo-autoupdate/autoupdate