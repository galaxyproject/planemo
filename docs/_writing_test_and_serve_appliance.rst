::

    $ planemo t
    ... Galaxy starts and runs the test ...
    All 1 test(s) executed passed.
    seqtk_seq[0]: passed

You can use the following command to open up the test results in your browser.

::

    $ firefox /opt/galaxy/tools/tool_test_output.html

Normally ``planemo`` requires an existing Galaxy instance to point at to run
the ``t`` (or ``test`` command) - but the virtual appliance has a
Galaxy instance preconfigured and registered with ``planemo``.
