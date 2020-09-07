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
  - ``--test``, which runs tests on all autoupdated tools
  - ``--update_test_data`` (if ``--test`` is also selected) which will update test data for failed tests

One of the most efficient ways to use it is by implementing a CI cron job which runs the command on an entire GitHub repository of tools.

Formatting tools
=============================

 ``autoupdate`` assumes that XML tools are formatted in a certain way - in accordance with the `IUC best practices`_

<https://planemo.readthedocs.io/en/latest/standards/docs/best_practices/tool_xml.html/>`__. In particular:

  - the tool ``version`` attribute must be set to ``@TOOL_VERSION@+galaxy0`` or ``@TOOL_VERSION@+galaxy@GALAXY_VERSION``
  - A token ``@TOOL_VERSION@`` should be created which corresponds to the version number of the main requirement.
  - Optionally, a token ``@GALAXY_VERSION@`` should be created, which should be an integer representing the number of times the XML wrapper has been updated since ``@TOOL_VERSION@`` was updated.


Implementing an autoupdate CI job
=================================

WIP

.. _IUC best practices: https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices/tool_xml.html
