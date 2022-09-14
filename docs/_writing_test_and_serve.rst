
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
``--galaxy_root`` to specify the root of a development Galaxy directory
to use (by default planemo will download and manage a instance just for planemo).

Finally, Planemo also allows us to combine the features of ``test`` and
``serve`` using the following command:

::

    $ planemo t --serve
    ...
    All 1 test(s) executed passed.
    ...

This command runs the specified tests, just like ``planemo test``, but
instead of shutting down the temporary Galaxy server afterwards, it
keeps it running, just like ``planemo serve``. You can open
http://127.0.0.1:9090 in a web browser and view the datasets generated
by the tests there. Each test is contained in a separate Galaxy history.