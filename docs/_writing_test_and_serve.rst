
::

    $ planemo t
    ...
    All 1 test(s) executed passed.
    seqtk_seq[0]: passed

.. include:: _writing_test_reports.rst

Now we can open Galaxy with the ``serve`` (or just ``s``).

::

    $ planemo s
    ...
    serving on http://127.0.0.1:9090

Open up http://127.0.0.1:9090 in a web browser to view your new tool.

Serve and test can be passed various command line arguments such as
``--galaxy_root`` to specify a Galaxy instance to use (by default
planemo will download and manage a instance just for planemo).

To inspect the histories and datasets produced by a test in Galaxy, combine
the two operations with ``test --serve``.

::

    $ planemo test --serve seqtk_seq.xml
    ...
    All 1 test(s) executed passed.
    Galaxy is serving the test histories at http://127.0.0.1:9090. Press Ctrl-C to stop.

Planemo writes the configured test reports before it waits, and keeps the same
Galaxy instance and its test histories available until it is interrupted. Use
``--host`` and ``--port`` to change the address. The option is supported by
Planemo-managed Galaxy engines, including the Gravity-backed
``installed_galaxy`` engine. It is not applicable to cwltool, Toil, or an
external Galaxy server whose lifecycle Planemo does not own.
