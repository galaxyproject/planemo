In addition to the in console display of test results as red (failing) or green
(passing), Planemo also creates an HTML report for the test results by default. Many
more test report options are available such ``--test_output_xunit`` which is useful
in certain `continuous integration`_ environments. See ``planemo test --help`` for
more options, as well as the ``test_reports`` `command
<http://planemo.readthedocs.io/en/latest/commands.html#test-reports-command>`__.

The machine-readable JSON report written by ``--test_output_json`` uses the
Planemo test report shape: a top-level ``version``, ``tests``, optional
``summary``, and optional ``exit_code``. ``test_reports`` validates this JSON
before rendering derived reports. ``merge_test_reports`` also validates each
input report and writes the same full report shape, including recalculated
``summary`` and ``exit_code`` fields. Older reports without ``summary`` are
still accepted and summarized during rendering.

.. _continuous integration: https://en.wikipedia.org/wiki/Continuous_integration
