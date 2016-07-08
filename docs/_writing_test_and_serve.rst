
::

    $ planemo t
    ...
    All 1 test(s) executed passed.
    seqtk_seq[0]: passed

Now we can open Galaxy with the ``serve`` (or just ``s``).

::

    $ planemo s
    ...
    serving on http://127.0.0.1:9090

Open up http://127.0.0.1:9090 in a web browser to view your new tool.

Serve and test can be passed various command line arguments such as
``--galaxy_root`` to specify a Galaxy instance to use (by default 
planemo will download and manage a instance just for planemo).
